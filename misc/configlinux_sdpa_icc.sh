#!/bin/sh

CC="icc"; export CC
CXX="icpc"; export CXX
FC="ifort"; export FC

./configure --prefix=$HOME/mplapack-work/MPLAPACK

