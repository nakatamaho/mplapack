#!/bin/bash

MPLIBS="_Float128 _Float64x dd double"

for _mplib in $MPLIBS; do
./Raxpy.${_mplib}_opt -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 10  >& log.Raxpy.${_mplib}_opt
./Raxpy.${_mplib} -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 3      >& log.Raxpy.${_mplib}
done


MPLIBS="mpfr gmp qd"

for _mplib in $MPLIBS; do
./Raxpy.${_mplib}_opt -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 10  >& log.Raxpy.${_mplib}_opt
./Raxpy.${_mplib} -NOCHECK -TOTALSTEPS 100 -STEP 3000 -LOOP 3       >& log.Raxpy.${_mplib}
done

gnuplot Raxpy1.plt > Raxpy1.eps
gnuplot Raxpy2.plt > Raxpy2.eps
