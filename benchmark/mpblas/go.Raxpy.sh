#!/bin/bash

####
MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
./Raxpy.${_mplib}_opt -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 10  >& log.Raxpy.${_mplib}_opt
./Raxpy.${_mplib} -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 3      >& log.Raxpy.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
./Raxpy.${_mplib}_opt -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 10  >& log.Raxpy.${_mplib}_opt
./Raxpy.${_mplib} -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 3       >& log.Raxpy.${_mplib}
done
####

####
if [ `uname` = "Linux" ]; then
    MODELNAME=`cat /proc/cpuinfo | grep 'model name' | uniq | awk '{for(i=4;i<=NF;i++) printf $i FS }'`
    SED=sed
elif [ `uname` = "Darwin" ]; then
    MODELNAME=`sysctl machdep.cpu.brand_string | awk '{for(i=2;i<=NF;i++) printf $i FS }'`
    SED=gsed
else
    MODELNAME="unknown"
fi

$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Raxpy1.plt
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Raxpy2.plt
####

gnuplot Raxpy1.plt > Raxpy1.eps
gnuplot Raxpy2.plt > Raxpy2.eps
