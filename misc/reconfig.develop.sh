#!/bin/sh

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77
ccache -M 20G
else
CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes --enable-debug=yes --enable-optimization=no

