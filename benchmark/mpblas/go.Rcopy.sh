#!/bin/bash

if [ `uname` = "Darwin" ]; then
    LDPATHPREFIX="DYLD_LIBRARY_PATH=%%PREFIX%%/lib"
else
    LDPATHPREFIX="LD_LIBRARY_PATH=%%PREFIX%%/lib:$LD_LIBRARY_PATH"
fi

####
MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rcopy.${_mplib}_opt -NOCHECK  >& log.Rcopy.${_mplib}_opt
env $LDPATHPREFIX ./Rcopy.${_mplib}     -NOCHECK  >& log.Rcopy.${_mplib}
done
####

####
MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
env $LDPATHPREFIX ./Rcopy.${_mplib}_opt -NOCHECK  >& log.Rcopy.${_mplib}_opt
env $LDPATHPREFIX ./Rcopy.${_mplib}     -NOCHECK  >& log.Rcopy.${_mplib}
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

$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rcopy1.plt
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rcopy2.plt
$SED -i -e "s/%%MODELNAME%%/$MODELNAME/g" Rcopy3.plt
####

gnuplot Rcopy1.plt > Rcopy1.pdf
gnuplot Rcopy2.plt > Rcopy2.pdf
gnuplot Rcopy3.plt > Rcopy3.pdf
