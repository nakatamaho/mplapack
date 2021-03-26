#!/bin/sh

CXX="/opt/local/bin/g++-mp-4.6" ; export CXX
CC="/opt/local/bin/gcc-mp-4.6" ; export CC
FC="/opt/local/bin/gfortran-mp-4.6"; export FC

./configure --prefix=$HOME/mplapack-work/MPLAPACK --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes --enable-_Float128=yes

