#!/bin/bash

Rpotrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3   >& log.Rpotrf.dd_cuda_total

gnuplot Rpotrf_cuda.plt > Rpotrf_cuda.eps
