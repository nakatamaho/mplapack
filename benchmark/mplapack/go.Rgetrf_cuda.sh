#!/bin/bash


if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX ./Rgetrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3   >& log.Rgetrf.dd_cuda_total

gnuplot Rgetrf_cuda.plt > Rgetrf_cuda.eps
