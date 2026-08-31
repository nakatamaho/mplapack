#!/bin/sh
# Build and run an external MPFR consumer using only pkg-config metadata.
# usage: check_mpfr_pkgconfig.sh <shared|static> [pkg-config-name] [source]
set -eu

if test "$#" -lt 1 || test "$#" -gt 3; then
    echo "usage: $0 <shared|static> [pkg-config-name] [source]" >&2
    exit 2
fi

mode=$1
pc_name=${2:-mplapack_mpfr}
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source=${3:-"$script_dir/mpfr_pkgconfig_consumer.cpp"}
pkg_config=${PKG_CONFIG:-pkg-config}
cxx=${CXX:-c++}

case "$mode" in
    shared) static_option= ;;
    static) static_option=--static ;;
    *) echo "invalid mode: $mode" >&2; exit 2 ;;
esac

if ! "$pkg_config" --exists "$pc_name"; then
    echo "FAIL: pkg-config package '$pc_name' is unavailable" >&2
    exit 1
fi

cflags=$("$pkg_config" --cflags $static_option "$pc_name")
libs=$("$pkg_config" --libs $static_option "$pc_name")
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-mpfr-pkgconfig.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
executable=$tmpdir/consumer

# CXX may intentionally contain a wrapper and its arguments, such as ccache.
# --static changes pkg-config's dependency closure, but does not force the
# primary package library itself to be selected from its archive.  Use the
# linker mode explicitly so this test really exercises static linkage.
if test "$mode" = static; then
    $cxx ${CXXFLAGS:-} "$source" $cflags -Wl,-Bstatic $libs -Wl,-Bdynamic -o "$executable"
else
    $cxx ${CXXFLAGS:-} "$source" $cflags $libs -o "$executable"
fi

runtime_dirs=
for flag in $libs; do
    case "$flag" in
        -L*)
            directory=${flag#-L}
            if test -n "$runtime_dirs"; then
                runtime_dirs=$runtime_dirs:$directory
            else
                runtime_dirs=$directory
            fi
            ;;
    esac
done

LD_LIBRARY_PATH="$runtime_dirs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
DYLD_LIBRARY_PATH="$runtime_dirs${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" \
PATH="$runtime_dirs${runtime_dirs:+:}$PATH" \
    "$executable"
echo "PASS: $mode pkg-config consumer for $pc_name"
