#!/bin/sh

CXX="g++-9" ; export CXX
CC="gcc-9" ; export CC
FC="gfortran-9"; export FC
F77="gfortran-9"; export F77

pushd mlapack/debug ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

#currently __float128 is not supported by gcc9 and 10. moreover dd builds but not give a correct answer.

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-longdouble=yes --enable-debug=yes





