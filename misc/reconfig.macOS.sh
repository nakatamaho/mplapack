#!/bin/bash

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache g++-mp-14" ; export CXX
CC="ccache gcc-mp-14" ; export CC
FC="gfortran-mp-14"; export FC
ccache -M 60G
else
CXX="g++-mp-14" ; export CXX
CC="gcc-mp-14" ; export CC
FC="gfortran-mp-14"; export FC
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
autoreconf --force --install

CONFIG_FLAGS=(
  "--prefix=$HOME/MPLAPACK"
  "--enable-gmp=yes"
  "--enable-mpfr=yes"
  "--enable-_Float128=yes"
  "--enable-qd=yes"
  "--enable-dd=yes"
  "--enable-double=yes"
  "--enable-test=yes"
  "--enable-benchmark=yes"
)

arch="$(uname -m)"
case "$arch" in
  x86_64|amd64)
    CONFIG_FLAGS+=("--enable-_Float64x=yes")
    ;;
  arm64|aarch64)
    CONFIG_FLAGS+=("--enable-_Float64x=no")
    ;;
  *)
    echo "Unknown arch: $arch (disabling _Float64x by default)" >&2
    CONFIG_FLAGS+=("--enable-_Float64x=no")
    ;;
esac

./configure "${CONFIG_FLAGS[@]}"
