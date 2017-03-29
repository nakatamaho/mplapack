#!/bin/sh

CXX="/opt/local/bin/g++-mp-4.6" ; export CXX
CC="/opt/local/bin/gcc-mp-4.6" ; export CC
FC="/opt/local/bin/gfortran-mp-4.6"; export FC

./configure --prefix=$HOME/mpack-work/MPACK --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes --enable-__float128=yes

