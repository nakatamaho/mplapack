#!/bin/bash

set -euo pipefail

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache clang++" ; export CXX
CC="ccache clang" ; export CC
FC="ccache flang-new"; export FC
ccache -M 80G
else
CXX="clang++" ; export CXX
CC="clang" ; export CC
FC="flang-new"; export FC
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
autoreconf --force --install

CONFIG_FLAGS=(
  "--prefix=$HOME/MPLAPACK"
  "--enable-gmp=yes"
  "--enable-mpfr=yes"
  "--enable-binary128=no"
  "--enable-qd=yes"
  "--enable-dd=no"
  "--enable-double=no"
  "--enable-test=yes"
  "--enable-benchmark=no"
  "--enable-optimization=no"
)

arch="$(uname -m)"
case "$arch" in
  x86_64|amd64)
    CONFIG_FLAGS+=("--enable-binary80=no")
    ;;
  arm64|aarch64)
    CONFIG_FLAGS+=("--enable-binary80=no")
    ;;
  *)
    echo "Unknown arch: $arch (disabling binary80 by default)" >&2
    CONFIG_FLAGS+=("--enable-binary80=no")
    ;;
esac

CFLAGS='-O1 -g -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls'
CXXFLAGS='-O1 -g -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls'
FCFLAGS='-O1 -g -fno-omit-frame-pointer'
LDFLAGS='-fsanitize=address'

CFLAGS="$CFLAGS" \
CXXFLAGS="$CXXFLAGS" \
FCFLAGS="$FCFLAGS" \
LDFLAGS="$LDFLAGS" \
./configure "${CONFIG_FLAGS[@]}"
