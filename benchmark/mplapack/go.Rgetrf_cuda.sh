#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX ./Rgetrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPM 5 -STEPN 5   >& log.Rgetrf.dd_cuda_total

gnuplot Rgetrf_cuda.plt > Rgetrf_cuda.pdf
