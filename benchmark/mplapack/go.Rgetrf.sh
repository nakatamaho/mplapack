#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX ./dgetrf_ref       -TOTALSTEPS 500 -STEPM 7 -STEPN 7 >& log.dgetrf.ref
env $LDPATHPREFIX ./dgetrf_openblas  -TOTALSTEPS 500 -STEPM 7 -STEPN 7 >& log.dgetrf.openblas

####
MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rgetrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPM 5 -STEPN 5    >& log.Rgetrf.${_mplib}_opt
env $LDPATHPREFIX ./Rgetrf.${_mplib}     -NOCHECK -TOTALSTEPS 200 -STEPM 5 -STEPN 5    >& log.Rgetrf.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rgetrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPM 5 -STEPN 5     >& log.Rgetrf.${_mplib}_opt
env $LDPATHPREFIX ./Rgetrf.${_mplib}     -NOCHECK -TOTALSTEPS 200 -STEPM 5 -STEPN 5     >& log.Rgetrf.${_mplib}
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

$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgetrf1.plt.in > Rgetrf1.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgetrf2.plt.in > Rgetrf2.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgetrf3.plt.in > Rgetrf3.plt
####

gnuplot Rgetrf1.plt > Rgetrf1.pdf
gnuplot Rgetrf2.plt > Rgetrf2.pdf
gnuplot Rgetrf3.plt > Rgetrf3.pdf
