#!/bin/sh

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache g++-10" ; export CXX
CC="ccache gcc-10" ; export CC
FC="gfortran-10"; export FC
ccache -M 80G
else
CXX="g++-10" ; export CXX
CC="gcc-10" ; export CC
FC="gfortran-10"; export FC
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-cuda=yes --enable-test=yes --enable-debug=yes --enable-optimization=no

