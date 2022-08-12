#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

####
MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rpotrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5   >& log.Rgemm.${_mplib}_opt
env $LDPATHPREFIX ./Rpotrf.${_mplib} -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3       >& log.Rgemm.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rpotrf.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5    >& log.Rgemm.${_mplib}_opt
env $LDPATHPREFIX ./Rpotrf.${_mplib} -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3        >& log.Rgemm.${_mplib}
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

$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rpotrf1.eps
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rpotrf2.eps
####

gnuplot Rpotrf1.plt > Rpotrf1.eps
gnuplot Rpotrf2.plt > Rpotrf2.eps
