#!/bin/bash

TIMEOUT=7200

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout $TIMEOUT ./dgemv_ref      -TOTALSTEPS 10000 >& log.dgemv.ref
env $LDPATHPREFIX timeout $TIMEOUT ./dgemv_openblas -TOTALSTEPS 10000 >& log.dgemv.openblas

####
MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout $TIMEOUT ./Rgemv.${_mplib}_opt -NOCHECK -TOTALSTEPS 10000 >& log.Rgemv.${_mplib}_opt
env $LDPATHPREFIX timeout $TIMEOUT ./Rgemv.${_mplib}     -NOCHECK -TOTALSTEPS 10000 >& log.Rgemv.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout $TIMEOUT ./Rgemv.${_mplib}_opt -NOCHECK  >& log.Rgemv.${_mplib}_opt
env $LDPATHPREFIX timeout $TIMEOUT ./Rgemv.${_mplib}     -NOCHECK  >& log.Rgemv.${_mplib}
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

$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgemv1.plt.in > Rgemv1.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgemv2.plt.in > Rgemv2.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Rgemv3.plt.in > Rgemv3.plt
####

gnuplot Rgemv1.plt > Rgemv1.pdf
gnuplot Rgemv2.plt > Rgemv2.pdf
gnuplot Rgemv3.plt > Rgemv3.pdf
