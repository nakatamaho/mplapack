#!/bin/sh

CXX="g++46" ; export CXX
CC="gcc46" ; export CC
FC="gfortran46"; export FC

CXXFLAGS="-O2 -I/usr/local/include" ; export CXXFLAGS
CPPFLAGS="-O2 -I/usr/local/include" ; export CPPFLAGS  
CFLAGS="-O2 -I/usr/local/include" ; export CFLAGS  
LDFLAGS="-L/usr/local/lib -Wl,-rpath=/usr/local/lib/gcc46"; export LDFLAGS
./configure --prefix=$HOME/mplapack-work-full/MPLAPACK --with-system-gmp=yes --with-system-mpfr=yes --with-system-mpc=yes --with-system-qd=yes --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes --enable-__float128=yes

