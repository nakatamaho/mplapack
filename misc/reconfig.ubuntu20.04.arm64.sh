#!/bin/sh

CXX="g++-9" ; export CXX
CC="gcc-9" ; export CC
FC="gfortran-9"; export FC
F77="gfortran-9"; export F77

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

#CXXFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CXXFLAGS
#CPPFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CPPFLAGS  
#CFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CFLAGS  

#currently __float128 is not supported by gcc9 and 10. dd builds but not give a correct answer.
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-__float128=no --enable-dd=no



