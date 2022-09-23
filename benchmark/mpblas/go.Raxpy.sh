#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

env $LDPATHPREFIX ./daxpy_ref                    >& log.daxpy.ref
env $LDPATHPREFIX ./daxpy_openblas               >& log.daxpy.openblas

####
MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Raxpy.${_mplib}_opt -NOCHECK  >& log.Raxpy.${_mplib}_opt
env $LDPATHPREFIX ./Raxpy.${_mplib} -NOCHECK      >& log.Raxpy.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Raxpy.${_mplib}_opt -NOCHECK  >& log.Raxpy.${_mplib}_opt
env $LDPATHPREFIX ./Raxpy.${_mplib} -NOCHECK      >& log.Raxpy.${_mplib}
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
