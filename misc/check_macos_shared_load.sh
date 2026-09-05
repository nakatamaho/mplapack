#!/bin/sh
# Build a disposable loader. Only the build-tree probe receives an RPATH.
# usage: check_macos_shared_load.sh <source> <build|installed> <dependency-dir> <library> ...
set -eu
test "$#" -ge 4 || exit 2
source=$1
mode=$2
dependency_dir=$3
shift 3
probe_dir=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-load.XXXXXX")
trap 'rm -rf "$probe_dir"' EXIT HUP INT TERM
case "$mode" in
    build) set -- "-Wl,-rpath,$dependency_dir" "$@" ;;
    installed) set -- "" "$@" ;;
    *) exit 2 ;;
esac
rpath=$1
shift
# CC may contain a launcher such as ccache. No backend library is linked.
if test -n "$rpath"; then
    ${CC:-cc} "$source" "$rpath" -o "$probe_dir/load"
else
    ${CC:-cc} "$source" -o "$probe_dir/load"
fi
# Ensure the installed probe cannot accidentally use build-tree overrides.
unset DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH DYLD_INSERT_LIBRARIES
for library in "$@"; do
    test -f "$library" || { echo "FAIL: missing library: $library" >&2; exit 1; }
    "$probe_dir/load" "$library"
done
