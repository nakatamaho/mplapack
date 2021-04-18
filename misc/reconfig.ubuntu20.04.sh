#!/bin/bash

USE_CCACHE=yes

if [ x"$USE_CCACHE" = x"yes" ] ; then
CXX="ccache g++-9" ; export CXX
CC="ccache gcc-9" ; export CC
FC="gfortran-9"; export FC
F77="gfortran-9"; export F77
ccache -M 20G
else
CXX="g++-9" ; export CXX
CC="gcc-9" ; export CC
FC="gfortran-9"; export FC
F77="gfortran-9"; export F77
fi

pushd mplapack/test ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
if [ `uname -m` = "x86_64" ]; then
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
else
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes
fi
