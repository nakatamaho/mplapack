#!/bin/bash
set -euo pipefail

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
    CXX="ccache clang++" ; export CXX
    CC="ccache clang" ; export CC
    FC="ccache gfortran-mp-15"; export FC
    ccache -M 80G
else
    CXX="clang++" ; export CXX
    CC="clang" ; export CC
    FC="gfortran-mp-15"; export FC
fi

# Generate Makefile.am and run autoreconf
pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
autoreconf --force --install

CONFIG_FLAGS=(
  "--prefix=$HOME/MPLAPACK"
  "--enable-gmp=no"
  "--enable-mpfr=yes"
  "--enable-qd=yes"
  "--enable-dd=no"
  "--enable-binary128=no"
  "--enable-binary80=no"
  "--enable-double=no"
  "--enable-test=yes"
  "--enable-benchmark=no"
  "--enable-optimization=no"
)

# Address Sanitizer (ASAN) and Debugging flags
ASAN_CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all"
ASAN_CXXFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all"
ASAN_LDFLAGS="-fsanitize=address,undefined"
FORTRAN_DEBUG_FLAGS="-O1 -g -fno-omit-frame-pointer"

# Run configure with exported environment variables
CFLAGS="${ASAN_CFLAGS}" \
CXXFLAGS="${ASAN_CXXFLAGS}" \
LDFLAGS="${ASAN_LDFLAGS}" \
FCFLAGS="${FORTRAN_DEBUG_FLAGS}" \
./configure "${CONFIG_FLAGS[@]}" "$@"
