#!/usr/bin/env bash

cd /home/docker/mplapack/mplapack/test/matgen

patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatm5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatmt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatms.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatms.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clarnd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Claror.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlarnd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatmt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlaror.cpp

cd /home/docker/mplapack/mplapack/test/lin/common
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkaa.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatb5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rtsqr01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrvls.cpp
