#!/bin/bash

TIMEOUT=7200

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout $TIMEOUT ./dsyrk_ref      -STEP 7 -TOTALSTEPS 714 -LOOPS 3 >& log.dsyrk.ref
env $LDPATHPREFIX timeout $TIMEOUT ./dsyrk_openblas -STEP 7 -TOTALSTEPS 714 -LOOPS 3 >& log.dsyrk.openblas

####
MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout $TIMEOUT ./Rsyrk.${_mplib}_opt -NOCHECK    >& log.Rsyrk.${_mplib}_opt
env $LDPATHPREFIX timeout $TIMEOUT ./Rsyrk.${_mplib}     -NOCHECK    >& log.Rsyrk.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout $TIMEOUT ./Rsyrk.${_mplib}_opt -NOCHECK    >& log.Rsyrk.${_mplib}_opt
env $LDPATHPREFIX timeout $TIMEOUT ./Rsyrk.${_mplib}     -NOCHECK    >& log.Rsyrk.${_mplib}
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

$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rsyrk1.plt.in > Rsyrk1.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rsyrk2.plt.in > Rsyrk2.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rsyrk3.plt.in > Rsyrk3.plt
####

gnuplot Rsyrk1.plt > Rsyrk1.pdf
gnuplot Rsyrk2.plt > Rsyrk2.pdf
gnuplot Rsyrk3.plt > Rsyrk3.pdf
