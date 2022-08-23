#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

####
MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rdot.${_mplib}_opt -NOCHECK  >& log.Rdot.${_mplib}_opt
env $LDPATHPREFIX ./Rdot.${_mplib}     -NOCHECK  >& log.Rdot.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rdot.${_mplib}_opt -NOCHECK  >& log.Rdot.${_mplib}_opt
env $LDPATHPREFIX ./Rdot.${_mplib}     -NOCHECK  >& log.Rdot.${_mplib}
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

$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rdot1.plt
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rdot2.plt
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rdot3.plt
####

gnuplot Rdot1.plt > Rdot1.pdf
gnuplot Rdot2.plt > Rdot2.pdf
gnuplot Rdot3.plt > Rdot3.pdf
