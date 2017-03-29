#!/bin/sh

CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
F77="ccache gfortran"; export F77

./configure --prefix=/work/MPACK
