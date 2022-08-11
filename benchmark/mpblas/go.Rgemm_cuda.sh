#!/bin/bash

./Rgemm.dd_cuda_kernel -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5   >& log.Rgemm.dd_cuda_kernel
./Rgemm.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3   >& log.Rgemm.dd_cuda_total

gnuplot Rgemm_cuda.plt > Rgemm_cuda.eps
