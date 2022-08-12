#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX ./Rpotrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEP 5     >& log.Rpotrf.dd_cuda_total

gnuplot Rpotrf_cuda.plt > Rpotrf_cuda.pdf
