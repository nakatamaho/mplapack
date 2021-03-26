#!/bin/sh

USE_CCACHE=yes

if [ x"$USE_CCACHE" = x"yes" ] ; then
CXX="ccache mingw64-g++" ; export CXX
CC="ccache mingw64-gcc" ; export CC
FC="mingw64-gfortran"; export FC
F77="mingw64-gfortran"; export F77
NM="mingw64-nm" ; export NM
RANLIB="mingw64-ranlib" ; export NM
AR="mingw64-ar" ; export NM
ccache -M 20G
else
CXX="mingw64-g++" ; export CXX
CC="mingw64-gcc" ; export CC
FC="mingw64-gfortran"; export FC
F77="mingw64-gfortran"; export F77
NM="mingw64-nm" ; export NM
RANLIB="mingw64-ranlib" ; export NM
AR="mingw64-ar" ; export NM
fi

pushd mplapack/debug ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

./configure --prefix=$HOME/MPLAPACK_MINGW64 --host=i386-pc-mingw64 --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-debug=yes


