#!/bin/sh

CC="icc"; export CC
CXX="icpc"; export CXX

./configure --prefix=$HOME/mplapack-icc-work/MPLAPACK
