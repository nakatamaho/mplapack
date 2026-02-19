#!/bin/bash

USE_CCACHE=yes

if [ x"$USE_CCACHE" = x"yes" ] ; then
CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
FC="gfortran"; export FC
ccache -M 80G
else
CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

autoreconf --force --install
if [ "`uname -m`" = "x86_64" ] || [ "`uname -m`" = "i686" ]; then
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-binary80=yes --enable-test=yes --enable-benchmark=yes
else
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes --enable-benchmark=yes
fi
