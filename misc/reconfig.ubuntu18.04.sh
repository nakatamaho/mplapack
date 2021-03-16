#!/bin/sh

CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77

pushd mplapack/debug ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-debug=yes


