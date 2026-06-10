#!/bin/bash

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache g++-mp-15" ; export CXX
CC="ccache gcc-mp-15" ; export CC
FC="gfortran-mp-15"; export FC
"$(cd "$(dirname "$0")" && pwd)/ccache_set_min_maxsize.sh" 80G
else
CXX="g++-mp-15" ; export CXX
CC="gcc-mp-15" ; export CC
FC="gfortran-mp-15"; export FC
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
autoreconf --force --install

CONFIG_FLAGS=(
  "--prefix=$HOME/MPLAPACK"
  "--enable-gmp=yes"
  "--enable-mpfr=yes"
  "--enable-binary128=yes"
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
