#!/bin/sh
# Build and run an external backend consumer using only pkg-config metadata.
# usage: check_backend_pkgconfig.sh <shared|static> <package> <source>
set -eu

if test "$#" -ne 3; then
    echo "usage: $0 <shared|static> <package> <source>" >&2
    exit 2
fi

mode=$1
pc_name=$2
source=$3
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

cflags=$($pkg_config --cflags $static_option "$pc_name")
libs=$($pkg_config --libs $static_option "$pc_name")
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-backend-pkgconfig.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
executable=$tmpdir/consumer

if test "$mode" = static; then
    # --static selects the dependency closure; -Bstatic also ensures that the
    # primary MPLAPACK archive is tested rather than its shared counterpart.
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
    "$executable"
echo "PASS: $mode pkg-config consumer for $pc_name"
