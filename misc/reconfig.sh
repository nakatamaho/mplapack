#!/bin/sh

CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
F77="ccache gfortran"; export F77

aclocal ; autoconf ; automake
autoreconf --force --install

#CXXFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CXXFLAGS
#CPPFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CPPFLAGS  
#CFLAGS="-O2 -Wall -Werror -pedantic-errors" ; export CFLAGS  

./configure --prefix=$HOME/MPLAPACK --enable-cypress

