#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

####
MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rgetrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5   >& log.Rgemm.${_mplib}_opt
env $LDPATHPREFIX ./Rgetrf.${_mplib} -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3       >& log.Rgemm.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rgetrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5    >& log.Rgemm.${_mplib}_opt
env $LDPATHPREFIX ./Rgetrf.${_mplib} -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3        >& log.Rgemm.${_mplib}
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

$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rgetrf1.plt
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rgetrf2.plt
####

gnuplot Rgetrf1.plt > Rgetrf1.pdf
gnuplot Rgetrf2.plt > Rgetrf2.pdf
