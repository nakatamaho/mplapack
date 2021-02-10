#!/bin/sh

export CXX=mingw32-g++
export CC=mingw32-gcc
export NM=mingw32-nm
export RANLIB=mingw32-ranlib
export AR=mingw32-ar

./configure --prefix=$HOME/mpack-mingw-work/MPACK --host=i386-pc-mingw32
