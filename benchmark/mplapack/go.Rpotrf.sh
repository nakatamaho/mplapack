#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout 3600 ./dpotrf_ref       -TOTALSTEPS 500 -STEP 7 >& log.dpotrf.ref
env $LDPATHPREFIX timeout 3600 ./dpotrf_openblas  -TOTALSTEPS 500 -STEP 7 >& log.dpotrf.openblas

####
MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout 3600 ./Rpotrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEP 5    >& log.Rpotrf.${_mplib}_opt
env $LDPATHPREFIX timeout 3600 ./Rpotrf.${_mplib}     -NOCHECK -TOTALSTEPS 200 -STEP 5    >& log.Rpotrf.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout 3600 ./Rpotrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEP 5    >& log.Rpotrf.${_mplib}_opt
env $LDPATHPREFIX timeout 3600 ./Rpotrf.${_mplib}     -NOCHECK -TOTALSTEPS 200 -STEP 5    >& log.Rpotrf.${_mplib}
done
####

####
if [ `uname` = "Linux" ]; then
    MODELNAME=`lscpu | grep 'Model name' | uniq | awk '{for(i=3;i<=NF;i++) printf $i FS }'`
    SED=sed
elif [ `uname` = "Darwin" ]; then
    MODELNAME=`sysctl machdep.cpu.brand_string | awk '{for(i=2;i<=NF;i++) printf $i FS }'`
    SED=gsed
else
    MODELNAME="unknown"
fi

$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rpotrf1.plt.in > Rpotrf1.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rpotrf2.plt.in > Rpotrf2.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rpotrf3.plt.in > Rpotrf3.plt
####

gnuplot Rpotrf1.plt > Rpotrf1.pdf
gnuplot Rpotrf2.plt > Rpotrf2.pdf
gnuplot Rpotrf3.plt > Rpotrf3.pdf
