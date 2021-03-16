#!/bin/bash
CXX="g++-mp-10" ; export CXX
CC="gcc-mp-10" ; export CC
FC="gfortran-mp-10"; export FC
F77="gfortran-mp-10"; export F77

pushd mplapack/debug ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-longdouble=yes --enable-debug=yes
