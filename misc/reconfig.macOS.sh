#!/bin/bash

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache g++-mp-10" ; export CXX
CC="ccache gcc-mp-10" ; export CC
FC="gfortran-mp-10"; export FC
INSTALL="/usr/bin/install"; export INSTALL
ccache -M 60G
else
CXX="g++-mp-10" ; export CXX
CC="gcc-mp-10" ; export CC
FC="gfortran-mp-10"; export FC
INSTALL="/usr/bin/install"; export INSTALL
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes --enable-benchmark=yes
