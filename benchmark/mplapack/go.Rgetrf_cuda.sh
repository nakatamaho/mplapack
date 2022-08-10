#!/bin/bash

Rgetrf.dd_cuda_total  -NOCHECK -TOTALSTEPS 700 -STEPK 5 -STEPM 5 -STEPN 5 -LOOP 3   >& log.Rgetrf.dd_cuda_total

gnuplot Rgetrf_cuda.plt > Rgetrf_cuda.eps
