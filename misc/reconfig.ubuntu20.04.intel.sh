#!/bin/bash

USE_CCACHE=yes

#https://gmplib.org/list-archives/gmp-bugs/2014-September/003526.html
if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache icpc" ; export CXX
CC="ccache icc" ; export CC
FC="ifort"; export FC
F77="ifort"; export F77
ccache -M 20G
else
CXX="icpc" ; export CXX
CC="icc" ; export CC
FC="ifort"; export FC
F77="ifort"; export F77
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
