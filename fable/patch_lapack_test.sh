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
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatb5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrqrt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrqrtp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatsp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatsy.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cpbt01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cpot01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cppt01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cpst01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatb5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrlqtp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrlqt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cunhr_col01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrvrf3.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatb4.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatb4.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkeq.cpp 
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clattp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clattr.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkeq.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlattr.cpp
#patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkaa.cpp
#patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkaa.cpp
