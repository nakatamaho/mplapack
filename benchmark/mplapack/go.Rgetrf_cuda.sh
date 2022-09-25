#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout 3600 ./Rgetrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPM 5 -STEPN 5   >& log.Rgetrf.dd_cuda_total

if [ `uname` = "Linux" ]; then
    MODELNAME=`nvidia-smi --query-gpu=name --format=csv | tail -1`
    SED=sed
else
    MODELNAME="unknown"
fi

$SED -e "s/%%GPU%%/$MODELNAME/g" Rgetrf_cuda.plt.in > Rgetrf_cuda.plt

gnuplot Rgetrf_cuda.plt > Rgetrf_cuda.pdf
