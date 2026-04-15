#!/usr/bin/env bash

cd /home/docker/mplapack/mplapack/test/matgen

patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clarnd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Claror.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatm5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatms.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatmt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlarnd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlaror.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatms.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatmt.cpp

cd /home/docker/mplapack/mplapack/test/lin/common
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkaa.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkeq.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkrfp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrlqt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrlqtp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrqrt.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cerrqrtp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatb4.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatb5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatsp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clatsy.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clattp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clattr.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clqt04.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clqt05.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cpbt01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cpot01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cppt01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cpst01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cqrt04.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cqrt05.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Ctsqr01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cunhr_col01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cunhr_col02.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkaa.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkeq.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkrfp.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrvrf3.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatb4.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatb5.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlattr.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rorhr_col01.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rorhr_col02.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rtsqr01.cpp

cd /home/docker/mplapack/mplapack/test/lin/
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Ctest.in
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rtest.in

cd /home/docker/mplapack/mplapack/test/eig/common
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkbd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkee.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cchkgg.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Ccsdts.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrges.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrges3.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrgsx.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrves.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrvsg.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrvsg2stg.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cget24.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cget52.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Clctsx.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cslect.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkbd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkec.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkee.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rchkgg.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rckcsd.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rcsdts.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrgsx.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrst2stg.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrves.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrvsg.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rget24.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlatb9.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rlctsx.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rslect.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Rdrgvx.cpp
patch  < ~/mplapack/fable/3.9.1/lapack/patch-Cdrgvx.cpp

cd /home/docker/mplapack/mplapack/test/eig
patch < ~/mplapack/fable/3.9.1/lapack/patch-Cbb.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-Ced.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-Red.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-Rgg.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-csd.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-gsv.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-nep.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-se2.in
patch < ~/mplapack/fable/3.9.1/lapack/patch-sep.in
