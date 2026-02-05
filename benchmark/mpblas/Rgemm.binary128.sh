#!/bin/bash
NUMOFTHREADS=12
export LD_LIBRARY_PATH=../../mpblas/optimized/_Float128/.libs:../../mpblas/reference/.libs/:/usr/local/gcc4.6.3/lib64/:$LD_LIBRARY_PATH 
OMP_NUM_THREADS=$NUMOFTHREADS ./Rgemm._Float128 -NN > log.Rgemm._Float128.$NUMOFTHREADS.NN
OMP_NUM_THREADS=$NUMOFTHREADS ./Rgemm._Float128 -NN -M0 50 -N0 100 -K0 200 > log.Rgemm._Float128.$NUMOFTHREADS.NN
OMP_NUM_THREADS=$NUMOFTHREADS ./Rgemm._Float128 -NT > log.Rgemm._Float128.$NUMOFTHREADS.NT
OMP_NUM_THREADS=$NUMOFTHREADS ./Rgemm._Float128 -NN > log.Rgemm._Float128.$NUMOFTHREADS.TN
OMP_NUM_THREADS=$NUMOFTHREADS ./Rgemm._Float128 -TT > log.Rgemm._Float128.$NUMOFTHREADS.TT
exit
