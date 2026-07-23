#!/bin/bash

set -euo pipefail

USE_CCACHE=yes
CLANG_VERSION=22

if [ x$USE_CCACHE = x"yes" ] ; then
    CXX="ccache clang++" ; export CXX
    CC="ccache clang" ; export CC
    FC="ccache gfortran-mp-15"; export FC
    "$(cd "$(dirname "$0")" && pwd)/ccache_set_min_maxsize.sh" 80G
else
    CXX="clang++" ; export CXX
    CC="clang" ; export CC
    FC="gfortran-mp-15"; export FC
fi

# Regenerate the autotools files
autoreconf --force --install

CONFIG_FLAGS=(
  "--prefix=$HOME/MPLAPACK"
  "--enable-gmp=yes"
  "--enable-mpfr=yes"
  "--enable-binary128=no"
  "--enable-qd=yes"
  "--enable-dd=yes"
  "--enable-double=yes"
  "--enable-test=yes"
  "--enable-benchmark=yes"
)

arch="$(uname -m)"
case "$arch" in
  x86_64|amd64)
    CONFIG_FLAGS+=("--enable-binary80=yes")
    ;;
  arm64|aarch64)
    CONFIG_FLAGS+=("--enable-binary80=no")
    ;;
  *)
    echo "Unknown arch: $arch (disabling binary80 by default)" >&2
    CONFIG_FLAGS+=("--enable-binary80=no")
    ;;
esac

./configure "${CONFIG_FLAGS[@]}"
