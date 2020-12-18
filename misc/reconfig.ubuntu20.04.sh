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

./configure --prefix=$HOME/MPACK

