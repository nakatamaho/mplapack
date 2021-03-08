#!/bin/sh

export CXX=mingw32-g++
export CC=mingw32-gcc
export FC=mingw32-gfortran
export NM=mingw32-nm
export RANLIB=mingw32-ranlib
export AR=mingw32-ar

./configure --prefix=$HOME/mplapack-mingw-work-full/MPLAPACK --enable-debug=yes --enable-mpfr=yes --enable-double=yes --host=i386-pc-mingw32

##--enable-__float128=yes --enable-qd=yes --enable-dd=yes

