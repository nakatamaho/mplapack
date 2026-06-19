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

tmp_meta_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_meta_dir"' EXIT

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

git_value() {
    local fallback="$1"; shift
    git -C "$PROJECT_ROOT" "$@" 2>/dev/null || printf '%s\n' "$fallback"
}

base_ref="$(git_value unknown rev-parse HEAD)"
base_ref_short="$(git_value unknown rev-parse --short=12 HEAD)"
branch="$(git_value unknown branch --show-current)"
if [[ -z "$branch" ]]; then
    branch="detached"
fi

tmp_status="$tmp_meta_dir/source-status.txt"
tmp_patch="$tmp_meta_dir/source.patch"
git -C "$PROJECT_ROOT" status --short > "$tmp_status" 2>/dev/null || true
git -C "$PROJECT_ROOT" diff --binary HEAD > "$tmp_patch" 2>/dev/null || : > "$tmp_patch"
patch_sha="$(sha256_file "$tmp_patch")"
if [[ ! -s "$tmp_patch" ]]; then
    patch_sha="none"
fi
if git -C "$PROJECT_ROOT" diff --quiet HEAD -- 2>/dev/null; then
    dirty="no"
else
    dirty="yes"
fi

make -C "$PROJECT_ROOT" dist
mkdir -p "$SNAPSHOT_DIR"
cp -f "$tmp_status" "$STATUS_FILE"
cp -f "$tmp_patch" "$PATCH_FILE"

tarball="$(find "$PROJECT_ROOT" -maxdepth 1 -type f \( -name 'mplapack-*.tar.gz' -o -name 'mplapack-*.tar.xz' -o -name 'mplapack-*.tar.bz2' \) -print0 | xargs -0 ls -t 2>/dev/null | head -n 1 || true)"
if [[ -z "$tarball" ]]; then
    echo "ERROR: make dist did not produce an mplapack tarball in $PROJECT_ROOT" >&2
    exit 1
fi

snapshot_tarball="$SNAPSHOT_DIR/$(basename "$tarball")"
cp -f "$tarball" "$snapshot_tarball"
tarball_sha="$(sha256_file "$snapshot_tarball")"
source_label="${base_ref_short}-${tarball_sha:0:12}"
if [[ "$dirty" == "yes" ]]; then
    source_label="${base_ref_short}-dirty-${tarball_sha:0:12}"
fi

{
    echo "MPLAPACK_SOURCE_MODE=dist"
    echo "MPLAPACK_SOURCE_BASE_REF=$base_ref"
    echo "MPLAPACK_SOURCE_BASE_REF_SHORT=$base_ref_short"
    echo "MPLAPACK_SOURCE_BRANCH=$branch"
    echo "MPLAPACK_SOURCE_DIRTY=$dirty"
    echo "MPLAPACK_SOURCE_PATCH_SHA256=$patch_sha"
    echo "MPLAPACK_SOURCE_TARBALL_SHA256=$tarball_sha"
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
    printf 'MPLAPACK_SOURCE_BASE_REF=%s\n' "$(shell_quote "$base_ref")"
} > "$INFO_FILE"

cat "$METADATA_FILE" >&2
echo "MPLAPACK_SOURCE_INFO_FILE=$INFO_FILE" >&2
