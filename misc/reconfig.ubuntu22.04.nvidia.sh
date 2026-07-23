#!/bin/bash

USE_CCACHE=yes

if [ x"$USE_CCACHE" = x"yes" ] ; then
CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
FC="ccache gfortran"; export FC
F77="ccache gfortran"; export F77
"$(cd "$(dirname "$0")" && pwd)/ccache_set_min_maxsize.sh" 80G
else
CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77
fi


autoreconf --force --install
if [ `uname -m` = "x86_64" ]; then
./configure --prefix=$HOME/MPLAPACK_NVIDIA --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-binary80=yes --enable-cuda=yes --with-opencl=/usr/local/cuda --enable-test=yes --enable-benchmark=yes
else
./configure --prefix=$HOME/MPLAPACK_NVIDIA --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-cuda=yes --with-opencl=/usr/local/cuda --enable-test=yes --enable-benchmark=yes
fi
