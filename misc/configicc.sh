#!/bin/sh

CC="icc"; export CC
CXX="icpc"; export CXX

./configure --prefix=$HOME/mpack-icc-work/MPACK
