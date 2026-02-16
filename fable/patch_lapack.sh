#!/usr/bin/env bash

cd /home/docker/mplapack/mplapack/reference

echo "Post-processing mplapack/reference/iMlaenv.cpp with fix_iMlaenv.py"
python $HOME/mplapack/fable/fix_iMlaenv.py iMlaenv.cpp
python $HOME/mplapack/fable/fix_iMlaenv.py iMparmq.cpp

#https://claude.ai/chat/ff991a40-a314-41e6-8699-d718148aec7e
#https://claude.ai/chat/238930a7-ded0-4d54-bc32-362acb7ff57e
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chbevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chbevx_2stage.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chbgvx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cheevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cheevx_2stage.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chpevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsbevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsbevx_2stage.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsbgvx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rspevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rstevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsyevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsyevx_2stage.cpp

#https://claude.ai/chat/f3c3ff73-215d-44ed-9dae-8e75808e06d9
#https://chatgpt.com/c/69929fe4-b980-83a7-8fd7-0988abca4785
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cgejsv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cgesvj.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rgejsv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rgesvj.cpp

patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cbbcsd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cggev.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cggev3.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cggevx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cgttrs.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cheequb.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chegvd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chgeqz.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Chpgvd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clacon.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clags2.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clahef_aa.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Claic1.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clargv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clarnv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clartg.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clasyf_aa.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clatbs.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clatps.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Clatrs.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cpotf2.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cpptrf.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Csyequb.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Csytri2x.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Ctfttp.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Ctgevc.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Ctpttf.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Ctrexc.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Ctrsyl.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Cunbdb.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rbbcsd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rbdsvdx.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlacon.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rladiv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlaln2.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlarnv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrb.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrk.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rlasy2.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rorbdb.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rspgvd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rstebz.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsyequb.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rsygvd.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rtfttp.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-Rtpttf.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-iMieeeck.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-iMlaenv.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-iMparam2stage.cpp
patch -p3 < ~/mplapack/fable/3.9.1/lapack/patch-iMparmq.cpp
