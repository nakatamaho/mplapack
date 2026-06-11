#!/bin/bash
set -euxo pipefail

echo '=== ccache stats (before) ==='
ccache -s || true

tarball=$(find /input -maxdepth 1 -name 'mplapack-*.tar.*' -print -quit)
test -n "$tarball"

rm -rf /work/mplapack-tarball
mkdir -p /work/mplapack-tarball
tar xf "$tarball" -C /work/mplapack-tarball
cd /work/mplapack-tarball/mplapack-*

ARCH=$(dpkg --print-architecture)
echo "Detected architecture: $ARCH"
COMMON_OPTS="--enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes --enable-benchmark=yes"
if [ "$ARCH" = "amd64" ] || [ "$ARCH" = "i386" ]; then
    CONFIGURE_OPTS="$COMMON_OPTS --enable-binary80=yes"
else
    CONFIGURE_OPTS="$COMMON_OPTS"
fi

echo "Configure options: $CONFIGURE_OPTS"
./configure $CONFIGURE_OPTS
make -j"$(nproc)"
make install

echo '=== ccache stats (after) ==='
ccache -s || true
