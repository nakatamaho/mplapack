#!/bin/bash
CXX="g++-mp-9" ; export CXX
CC="gcc-mp-9" ; export CC
FC="gfortran-mp-9"; export FC
F77="gfortran-mp-9"; export F77

pushd mplapack/debug ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --with-system-gmp=/opt/local --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-debug=yes
