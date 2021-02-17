#!/bin/bash

CXX="g++-9" ; export CXX
CC="gcc-9" ; export CC
FC="gfortran-9"; export FC
F77="gfortran-9"; export F77

pushd mlapack/debug ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-__float128=yes --enable-debug=yes


