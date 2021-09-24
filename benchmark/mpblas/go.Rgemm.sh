#!/bin/bash

MPLIBS="_Float64x dd double _Float128"

for _mplib in $MPLIBS; do
./Rgemm.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5   >& log.Rgemm.${_mplib}_opt
./Rgemm.${_mplib} -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3       >& log.Rgemm.${_mplib}
done


MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
./Rgemm.${_mplib}_opt -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5    >& log.Rgemm.${_mplib}_opt
./Rgemm.${_mplib} -NOCHECK -TOTALSTEPS 200 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3        >& log.Rgemm.${_mplib}
done

gnuplot Rgemm1.plt > Rgemm1.eps
gnuplot Rgemm2.plt > Rgemm2.eps
