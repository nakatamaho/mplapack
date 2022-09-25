#!/bin/bash

TIMEOUT=7200

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout $TIMEOUT ./Rpotrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEP 5     >& log.Rpotrf.dd_cuda_total

if [ `uname` = "Linux" ]; then
    MODELNAME=`nvidia-smi --query-gpu=name --format=csv | tail -1`
    SED=sed
else
    MODELNAME="unknown"
fi

$SED -e "s/%%GPU%%/$MODELNAME/g" Rpotrf_cuda.plt.in > Rpotrf_cuda.plt

gnuplot Rpotrf_cuda.plt > Rpotrf_cuda.pdf
