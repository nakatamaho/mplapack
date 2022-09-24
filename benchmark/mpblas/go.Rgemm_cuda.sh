#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

####

env $LDPATHPREFIX ./Rgemm.dd_cuda_kernel -NOCHECK -TOTALSTEPS 720 -STEPK 7 -STEPM 7 -STEPN 7 -LOOP 3 >& log.Rgemm.dd_cuda_kernel
env $LDPATHPREFIX ./Rgemm.dd_cuda_total  -NOCHECK -TOTALSTEPS 720 -STEPK 7 -STEPM 7 -STEPN 7 -LOOP 3 >& log.Rgemm.dd_cuda_total

if [ `uname` = "Linux" ]; then
    MODELNAME=`nvidia-smi --query-gpu=name --format=csv | tail -1`
    SED=sed
else
    MODELNAME="unknown"
fi

$SED -e "s/%%GPU%%/$MODELNAME/g" Rgemm_cuda.plt.in > Rgemm_cuda.plt

gnuplot Rgemm_cuda.plt > Rgemm_cuda.pdf
