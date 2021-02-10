#!/bin/sh

CC="icc"; export CC
CXX="icpc"; export CXX
FC="ifort"; export FC

./configure --prefix=$HOME/mpack-icc-work-full/MPACK --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes


