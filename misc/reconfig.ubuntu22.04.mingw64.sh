#!/bin/sh

USE_CCACHE=yes

if [ x"$USE_CCACHE" = x"yes" ] ; then
CXX="ccache x86_64-w64-mingw32-g++" ; export CXX
CC="ccache x86_64-w64-mingw32-gcc" ; export CC
FC="x86_64-w64-mingw32-gfortran"; export FC
NM="x86_64-w64-mingw32-nm" ; export NM
RANLIB="x86_64-w64-mingw32-ranlib" ; export RANLIB
AR="x86_64-w64-mingw32-ar" ; export AR
ccache -M 80G
else
CXX="x86_64-w64-mingw32-g++" ; export CXX
CC="x86_64-w64-mingw32-gcc" ; export CC
FC="x86_64-w64-mingw32-gfortran"; export FC
NM="x86_64-w64-mingw32-nm" ; export NM
RANLIB="x86_64-w64-mingw32-ranlib" ; export RANLIB
AR="x86_64-w64-mingw32-ar" ; export AR
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK_MINGW --host=x86_64-w64-mingw32 --target=x86_64-w64-mingw32 --enable-gmp=yes --enable-mpfr=yes --enable-double=yes --enable-_Float64x=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-test=yes --enable-benchmark=yes



