#!/bin/bash

MPLIBS="_Float128 _Float64x dd double gmp mpfr qd"

for _mplib in $MPLIBS; do
./Raxpy.${_mplib}_opt -NOCHECK -STEP 400 -LOOP 10  >& log.Raxpy.${_mplib}_opt
./Raxpy.${_mplib} -NOCHECK -STEP 400 -LOOP 10      >& log.Raxpy.${_mplib}
done

