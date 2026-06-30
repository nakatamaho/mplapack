#!/bin/bash

mplapack_dist_cache_sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "ERROR: no SHA256 tool found (sha256sum or shasum)" >&2
        return 1
    fi
}

mplapack_dist_cache_sha256_stdin() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 | awk '{print $1}'
    else
        echo "ERROR: no SHA256 tool found (sha256sum or shasum)" >&2
        return 1
    fi
}

mplapack_dist_cache_key() {
    local project_root="$1"
    local base_ref patch_sha

    if ! git -C "$project_root" rev-parse --verify -q HEAD >/dev/null 2>&1; then
        return 1
    fi

    base_ref="$(git -C "$project_root" rev-parse HEAD)"
    if git -C "$project_root" diff --quiet HEAD -- 2>/dev/null; then
        printf 'git-%s\n' "$base_ref"
    else
        patch_sha="$(git -C "$project_root" diff --binary HEAD | mplapack_dist_cache_sha256_stdin)"
        printf 'git-%s-dirty-%s\n' "$base_ref" "${patch_sha:0:16}"
    fi
}

mplapack_dist_cache_find_tarball() {
    local cache_dir="$1"
    local cache_key="$2"
    local entry="$cache_dir/$cache_key"
    local candidate newest=""

    [[ -d "$entry" ]] || return 1
    while IFS= read -r -d '' candidate; do
        if [[ -z "$newest" || "$candidate" -nt "$newest" ]]; then
            newest="$candidate"
        fi
    done < <(find "$entry" -maxdepth 1 -type f \
        \( -name 'mplapack-*.tar.gz' -o -name 'mplapack-*.tar.xz' -o -name 'mplapack-*.tar.bz2' \) \
        -print0)
    [[ -n "$newest" ]] || return 1
    printf '%s\n' "$newest"
}

mplapack_dist_cache_reuse() {
    local cache_dir="$1"
    local cache_key="$2"
    local dest_dir="$3"
    local cached_tarball dest_tarball metadata_file expected_sha actual_sha

    if [[ "${MPLAPACK_DIST_CACHE_REUSE:-yes}" != "yes" ]]; then
        return 1
    fi

    cached_tarball="$(mplapack_dist_cache_find_tarball "$cache_dir" "$cache_key" || true)"
    [[ -n "$cached_tarball" && -f "$cached_tarball" ]] || return 1

    metadata_file="$(dirname "$cached_tarball")/dist-cache.env"
    [[ -f "$metadata_file" ]] || return 1
    expected_sha="$(awk -F= '$1 == "MPLAPACK_DIST_CACHE_TARBALL_SHA256" { print $2; exit }' "$metadata_file")"
    [[ -n "$expected_sha" ]] || return 1
    actual_sha="$(mplapack_dist_cache_sha256_file "$cached_tarball")"
    [[ "$actual_sha" == "$expected_sha" ]] || return 1

    mkdir -p "$dest_dir"
    dest_tarball="$dest_dir/$(basename "$cached_tarball")"
    cp -p "$cached_tarball" "$dest_tarball"
    if [[ -f "${cached_tarball}.sha256sum" ]]; then
        cp -p "${cached_tarball}.sha256sum" "${dest_tarball}.sha256sum"
    fi
    if [[ -f "${cached_tarball}.md5sum" ]]; then
        cp -p "${cached_tarball}.md5sum" "${dest_tarball}.md5sum"
    fi
    printf '%s\n' "$dest_tarball"
}

mplapack_dist_cache_store() {
    local cache_dir="$1"
    local cache_key="$2"
    local tarball="$3"
    local entry="$cache_dir/$cache_key"
    local dest_tarball tarball_sha

    [[ -n "$cache_key" && -f "$tarball" ]] || return 0

    mkdir -p "$entry"
    dest_tarball="$entry/$(basename "$tarball")"
    cp -p "$tarball" "$dest_tarball"
    if [[ -f "${tarball}.sha256sum" ]]; then
        cp -p "${tarball}.sha256sum" "${dest_tarball}.sha256sum"
    fi
    if [[ -f "${tarball}.md5sum" ]]; then
        cp -p "${tarball}.md5sum" "${dest_tarball}.md5sum"
    fi
    tarball_sha="$(mplapack_dist_cache_sha256_file "$dest_tarball")"
    {
        printf 'MPLAPACK_DIST_CACHE_KEY=%s\n' "$cache_key"
        printf 'MPLAPACK_DIST_CACHE_TARBALL=%s\n' "$(basename "$dest_tarball")"
        printf 'MPLAPACK_DIST_CACHE_TARBALL_SHA256=%s\n' "$tarball_sha"
        printf 'MPLAPACK_DIST_CACHE_CREATED_AT=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "$entry/dist-cache.env"
}
