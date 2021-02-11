#!/bin/sh

CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

#CXXFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CXXFLAGS
#CPPFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CPPFLAGS  
#CFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CFLAGS  

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-__float128=yes


