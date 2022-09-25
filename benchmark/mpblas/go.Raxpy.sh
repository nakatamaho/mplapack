#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX timeout 3600 ./daxpy_ref      -STEP 7 -TOTALSTEPS 42857 -LOOPS 3 >& log.daxpy.ref
env $LDPATHPREFIX timeout 3600 ./daxpy_openblas -STEP 7 -TOTALSTEPS 42857 -LOOPS 3 >& log.daxpy.openblas

####
MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout 3600 ./Raxpy.${_mplib}_opt -NOCHECK  >& log.Raxpy.${_mplib}_opt
env $LDPATHPREFIX timeout 3600 ./Raxpy.${_mplib} -NOCHECK      >& log.Raxpy.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX timeout 3600 ./Raxpy.${_mplib}_opt -NOCHECK  >& log.Raxpy.${_mplib}_opt
env $LDPATHPREFIX timeout 3600 ./Raxpy.${_mplib} -NOCHECK      >& log.Raxpy.${_mplib}
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

$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Raxpy1.plt.in > Raxpy1.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Raxpy2.plt.in > Raxpy2.plt
$SED -e "s/%%MODELNAME%%/$MODELNAME/g" Raxpy3.plt.in > Raxpy3.plt
####

gnuplot Raxpy1.plt > Raxpy1.pdf
gnuplot Raxpy2.plt > Raxpy2.pdf
gnuplot Raxpy3.plt > Raxpy3.pdf
