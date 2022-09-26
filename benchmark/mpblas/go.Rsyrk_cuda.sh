#!/bin/bash

TIMEOUT=7200

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout $TIMEOUT ./Rsyrk.dd_cuda_kernel -NOCHECK -TOTALSTEPS 720 -STEPK 7 -STEPN 7 -LOOP 3   >& log.Rsyrk.dd_cuda_kernel
env $LDPATHPREFIX timeout $TIMEOUT ./Rsyrk.dd_cuda_total  -NOCHECK -TOTALSTEPS 720 -STEPK 7 -STEPN 7 -LOOP 3   >& log.Rsyrk.dd_cuda_total

if [ `uname` = "Linux" ]; then
    MODELNAME=`nvidia-smi --query-gpu=name --format=csv | tail -1`
    SED=sed
else
    MODELNAME="unknown"
    SED=sed
fi

$SED -e "s/%%GPU%%/$MODELNAME/g" Rsyrk_cuda.plt.in > Rsyrk_cuda.plt

gnuplot Rsyrk_cuda.plt > Rsyrk_cuda.pdf
