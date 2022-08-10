#!/bin/bash

Rsyrk.dd_cuda_kernel -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 5   >& log.Rsyrk.dd_cuda_kernel
Rsyrk.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3   >& log.Rsyrk.dd_cuda_total

gnuplot Rsyrk_cuda.plt > Rsyrk_cuda.eps
