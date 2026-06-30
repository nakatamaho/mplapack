#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs/$(date +%Y%m%d_%H%M%S)}"
SNAPSHOT_DIR="${MPLAPACK_SOURCE_SNAPSHOT_DIR:-$LOGDIR/source}"
INFO_FILE="${MPLAPACK_SOURCE_INFO_FILE:-$SNAPSHOT_DIR/source.env}"
METADATA_FILE="$SNAPSHOT_DIR/source-metadata.txt"
PATCH_FILE="$SNAPSHOT_DIR/source.patch"
STATUS_FILE="$SNAPSHOT_DIR/source-status.txt"
SNAPSHOT_REUSE="${MPLAPACK_SOURCE_SNAPSHOT_REUSE:-yes}"
DIST_CACHE_DIR="${MPLAPACK_DIST_CACHE_DIR:-$SCRIPT_DIR/dist-cache}"

# shellcheck disable=SC1091
source "$SCRIPT_DIR/dist-cache.sh"

tmp_meta_dir="$(mktemp -d)"
DIST_ARTIFACT_DIR="$tmp_meta_dir/top-level-dist-artifacts"
DIST_RESTORE_NEEDED=no

find_top_level_dist_artifacts() {
    find "$PROJECT_ROOT" -maxdepth 1 -type f \
        \( -name 'mplapack-*.tar.gz' -o -name 'mplapack-*.tar.xz' -o -name 'mplapack-*.tar.bz2' \
        -o -name 'mplapack-*.tar.gz.sha256sum' -o -name 'mplapack-*.tar.xz.sha256sum' -o -name 'mplapack-*.tar.bz2.sha256sum' \
        -o -name 'mplapack-*.tar.gz.md5sum' -o -name 'mplapack-*.tar.xz.md5sum' -o -name 'mplapack-*.tar.bz2.md5sum' \) \
        -print0
}

preserve_top_level_dist_artifacts() {
    mkdir -p "$DIST_ARTIFACT_DIR"
    while IFS= read -r -d '' artifact; do
        cp -p "$artifact" "$DIST_ARTIFACT_DIR/"
    done < <(find_top_level_dist_artifacts)
    DIST_RESTORE_NEEDED=yes
}

restore_top_level_dist_artifacts() {
    if [[ "$DIST_RESTORE_NEEDED" != "yes" ]]; then
        return 0
    fi

    find_top_level_dist_artifacts | xargs -0 rm -f
    if compgen -G "$DIST_ARTIFACT_DIR/*" >/dev/null; then
        cp -p "$DIST_ARTIFACT_DIR"/* "$PROJECT_ROOT/"
    fi
    DIST_RESTORE_NEEDED=no
}

cleanup() {
    restore_top_level_dist_artifacts
    rm -rf "$tmp_meta_dir"
}
trap cleanup EXIT

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "ERROR: no SHA256 tool found (sha256sum or shasum)" >&2
        return 1
    fi
}

shell_quote() {
    printf '%q' "$1"
}

reuse_existing_snapshot() {
    if [[ "$SNAPSHOT_REUSE" != "yes" || ! -f "$INFO_FILE" ]]; then
        return 1
    fi

    # shellcheck disable=SC1090
    source "$INFO_FILE"
    if [[ "${MPLAPACK_SOURCE_MODE:-}" != "dist" || -z "${MPLAPACK_SOURCE_TARBALL:-}" || ! -f "$MPLAPACK_SOURCE_TARBALL" ]]; then
        echo "ERROR: existing source snapshot is incomplete: $INFO_FILE" >&2
        echo "       Remove $SNAPSHOT_DIR or set MPLAPACK_SOURCE_SNAPSHOT_REUSE=no to regenerate." >&2
        exit 1
    fi

    if [[ -n "${MPLAPACK_SOURCE_METADATA:-}" && -f "$MPLAPACK_SOURCE_METADATA" ]]; then
        cat "$MPLAPACK_SOURCE_METADATA" >&2
    fi
    echo "Reusing existing MPLAPACK_SOURCE_INFO_FILE=$INFO_FILE" >&2
    echo "MPLAPACK_SOURCE_INFO_FILE=$INFO_FILE" >&2
    exit 0
}

git_value() {
    local fallback="$1"; shift
    git -C "$PROJECT_ROOT" "$@" 2>/dev/null || printf '%s\n' "$fallback"
}

reuse_existing_snapshot || true

base_ref="$(git_value unknown rev-parse HEAD)"
base_ref_short="$(git_value unknown rev-parse --short=12 HEAD)"
branch="$(git_value unknown branch --show-current)"
if [[ -z "$branch" ]]; then
    branch="detached"
fi
if git -C "$PROJECT_ROOT" rev-parse --verify -q HEAD >/dev/null 2>&1; then
    git_head=yes
else
    git_head=no
fi

tmp_status="$tmp_meta_dir/source-status.txt"
tmp_patch="$tmp_meta_dir/source.patch"
git -C "$PROJECT_ROOT" status --short > "$tmp_status" 2>/dev/null || true
if [[ "$git_head" == "yes" ]]; then
    git -C "$PROJECT_ROOT" diff --binary HEAD > "$tmp_patch" 2>/dev/null || : > "$tmp_patch"
else
    : > "$tmp_patch"
fi
patch_sha="$(sha256_file "$tmp_patch")"
if [[ ! -s "$tmp_patch" ]]; then
    patch_sha="none"
fi
if [[ "$git_head" == "yes" ]] && git -C "$PROJECT_ROOT" diff --quiet HEAD -- 2>/dev/null; then
    dirty="no"
else
    dirty="yes"
fi

dist_cache_key="$(mplapack_dist_cache_key "$PROJECT_ROOT" 2>/dev/null || true)"
mkdir -p "$SNAPSHOT_DIR"
cp -f "$tmp_status" "$STATUS_FILE"
cp -f "$tmp_patch" "$PATCH_FILE"

snapshot_tarball=""
if [[ -n "$dist_cache_key" ]]; then
    snapshot_tarball="$(mplapack_dist_cache_reuse "$DIST_CACHE_DIR" "$dist_cache_key" "$SNAPSHOT_DIR" || true)"
    if [[ -n "$snapshot_tarball" ]]; then
        echo "Reusing cached dist tarball: $snapshot_tarball" >&2
        echo "Dist cache key: $dist_cache_key" >&2
    fi
fi

if [[ -z "$snapshot_tarball" ]]; then
    preserve_top_level_dist_artifacts
    make -C "$PROJECT_ROOT" dist

    tarball="$(find "$PROJECT_ROOT" -maxdepth 1 -type f \( -name 'mplapack-*.tar.gz' -o -name 'mplapack-*.tar.xz' -o -name 'mplapack-*.tar.bz2' \) -print0 | xargs -0 ls -t 2>/dev/null | head -n 1 || true)"
    if [[ -z "$tarball" ]]; then
        echo "ERROR: make dist did not produce an mplapack tarball in $PROJECT_ROOT" >&2
        exit 1
    fi

    snapshot_tarball="$SNAPSHOT_DIR/$(basename "$tarball")"
    cp -f "$tarball" "$snapshot_tarball"
    restore_top_level_dist_artifacts
    if [[ -n "$dist_cache_key" ]]; then
        mplapack_dist_cache_store "$DIST_CACHE_DIR" "$dist_cache_key" "$snapshot_tarball"
        echo "Stored dist tarball cache: $DIST_CACHE_DIR/$dist_cache_key" >&2
    fi
fi

tarball_sha="$(sha256_file "$snapshot_tarball")"
if [[ "$git_head" == "yes" ]]; then
    source_label="${base_ref_short}-${tarball_sha:0:12}"
    if [[ "$dirty" == "yes" ]]; then
        source_label="${base_ref_short}-dirty-${tarball_sha:0:12}"
    fi
else
    source_label="nogit-${tarball_sha:0:12}"
fi

{
    echo "MPLAPACK_SOURCE_MODE=dist"
    echo "MPLAPACK_SOURCE_BASE_REF=$base_ref"
    echo "MPLAPACK_SOURCE_BASE_REF_SHORT=$base_ref_short"
    echo "MPLAPACK_SOURCE_BRANCH=$branch"
    echo "MPLAPACK_SOURCE_GIT_HEAD=$git_head"
    echo "MPLAPACK_SOURCE_DIRTY=$dirty"
    echo "MPLAPACK_SOURCE_PATCH_SHA256=$patch_sha"
    echo "MPLAPACK_SOURCE_TARBALL_SHA256=$tarball_sha"
    echo "MPLAPACK_DIST_CACHE_KEY=${dist_cache_key:-}"
    echo "MPLAPACK_SOURCE_TARBALL=$(basename "$snapshot_tarball")"
    echo "MPLAPACK_SOURCE_LABEL=$source_label"
    echo "MPLAPACK_SOURCE_CREATED_AT=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$METADATA_FILE"

{
    printf 'MPLAPACK_SOURCE_MODE=dist\n'
    printf 'MPLAPACK_SOURCE_TARBALL=%s\n' "$(shell_quote "$snapshot_tarball")"
    printf 'MPLAPACK_SOURCE_METADATA=%s\n' "$(shell_quote "$METADATA_FILE")"
    printf 'MPLAPACK_SOURCE_PATCH=%s\n' "$(shell_quote "$PATCH_FILE")"
    printf 'MPLAPACK_SOURCE_STATUS=%s\n' "$(shell_quote "$STATUS_FILE")"
    printf 'MPLAPACK_SOURCE_LABEL=%s\n' "$(shell_quote "$source_label")"
    printf 'MPLAPACK_SOURCE_TARBALL_SHA256=%s\n' "$(shell_quote "$tarball_sha")"
    printf 'MPLAPACK_DIST_CACHE_KEY=%s\n' "$(shell_quote "${dist_cache_key:-}")"
    printf 'MPLAPACK_SOURCE_BASE_REF=%s\n' "$(shell_quote "$base_ref")"
    printf 'MPLAPACK_SOURCE_GIT_HEAD=%s\n' "$(shell_quote "$git_head")"
    printf 'MPLAPACK_SOURCE_DIRTY=%s\n' "$(shell_quote "$dirty")"
} > "$INFO_FILE"

cat "$METADATA_FILE" >&2
echo "MPLAPACK_SOURCE_INFO_FILE=$INFO_FILE" >&2
