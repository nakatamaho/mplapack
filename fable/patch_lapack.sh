#!/usr/bin/env bash

cd /home/docker/mplapack/mplapack/reference

echo "Post-processing mplapack/reference/iMlaenv.cpp with fix_iMlaenv.py"
python $HOME/mplapack/fable/fix_iMlaenv.py iMlaenv.cpp
python $HOME/mplapack/fable/fix_iMlaenv.py iMparmq.cpp

patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cbbcsd.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cpptrf.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cpotf2.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clarnv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Ctrexc.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clartg.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Rbbcsd.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrb.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlarnv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rstebz.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrk.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Claic1.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Csyequb.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cgesvj.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Rgesvj.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlacon.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clacon.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlasy2.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cgejsv.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Rgejsv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlaln2.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-iMieeeck.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rladiv.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cgttrs.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cggev3.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cheequb.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cggev.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cggevx.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clags2.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clartg.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clargv.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Clasyf_aa.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Clatbs.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-iMlaenv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-iMparmq.cpp
