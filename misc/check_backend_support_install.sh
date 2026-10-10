#!/bin/sh
# Exercise first-time support-library installation and libtool relinking.
# Usage: check_backend_support_install.sh <configured-top-builddir>
# Build the enabled primary and support libraries before running this check.
set -eu
if test "$#" -ne 1; then
    echo "usage: $0 <configured-top-builddir>" >&2
    exit 2
fi
builddir=$(cd "$1" && pwd)
stage=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-support-install.XXXXXX")
trap 'rm -rf "$stage"' EXIT HUP INT TERM
make_cmd=${MAKE:-make}

install_libraries() {
    directory=$1
    shift
    # Test installation of already built artifacts, without recompiling them
    # because configure-generated headers have newer timestamps.
    for library in "$directory"/*.la; do
        test -f "$library" || continue
        set -- "$@" -o "${library##*/}"
    done
    "$make_cmd" -C "$directory" "$@" DESTDIR="$stage" install-libLTLIBRARIES
}

# Stage prerequisites first, then exercise Automake's actual installation
# order. DESTDIR isolates the test from the configured installation prefix.
install_libraries "$builddir/mplapack/reference"
for directory in "$builddir"/mpblas/optimized/*; do
    backend=${directory##*/}
    if test -f "$directory/libmplapack_${backend}_opt.la"; then
        install_libraries "$directory"
    fi
done
for family in matgen lin eig; do
    install_libraries "$builddir/mplapack/test/$family"
done
echo "PASS: support libraries installed into a fresh DESTDIR"
