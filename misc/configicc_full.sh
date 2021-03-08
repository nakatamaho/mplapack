#!/bin/sh

CC="icc"; export CC
CXX="icpc"; export CXX
FC="ifort"; export FC

./configure --prefix=$HOME/mplapack-icc-work-full/MPLAPACK --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes


