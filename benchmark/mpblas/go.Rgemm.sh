#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout 3600 ./dgemm_ref      -STEP 7 -TOTALSTEPS 714 -LOOPS 3 >& log.dgemm.ref
env $LDPATHPREFIX timeout 3600 ./dgemm_openblas -STEP 7 -TOTALSTEPS 714 -LOOPS 3 >& log.dgemm.openblas

####
MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout 3600 ./Rgemm.${_mplib}_opt -NOCHECK    >& log.Rgemm.${_mplib}_opt
env $LDPATHPREFIX timeout 3600 ./Rgemm.${_mplib}     -NOCHECK    >& log.Rgemm.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout 3600 ./Rgemm.${_mplib}_opt -NOCHECK     >& log.Rgemm.${_mplib}_opt
env $LDPATHPREFIX timeout 3600 ./Rgemm.${_mplib}     -NOCHECK     >& log.Rgemm.${_mplib}
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

$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgemm1.plt.in > Rgemm1.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgemm2.plt.in > Rgemm2.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgemm3.plt.in > Rgemm3.plt
####

gnuplot Rgemm1.plt > Rgemm1.pdf
gnuplot Rgemm2.plt > Rgemm2.pdf
gnuplot Rgemm3.plt > Rgemm3.pdf
