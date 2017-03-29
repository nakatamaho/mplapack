#!/bin/bash

OMP_NUM_THREADS=64 ./Rgemm.mpfr > log.Rgemm.mpfr.64
OMP_NUM_THREADS=60 ./Rgemm.mpfr > log.Rgemm.mpfr.60
OMP_NUM_THREADS=48 ./Rgemm.mpfr > log.Rgemm.mpfr.48
OMP_NUM_THREADS=32 ./Rgemm.mpfr > log.Rgemm.mpfr.32
OMP_NUM_THREADS=30 ./Rgemm.mpfr > log.Rgemm.mpfr.30
OMP_NUM_THREADS=24 ./Rgemm.mpfr > log.Rgemm.mpfr.24
OMP_NUM_THREADS=16 ./Rgemm.mpfr > log.Rgemm.mpfr.16
OMP_NUM_THREADS=8 ./Rgemm.mpfr > log.Rgemm.mpfr.8
OMP_NUM_THREADS=4 ./Rgemm.mpfr > log.Rgemm.mpfr.4
OMP_NUM_THREADS=2 ./Rgemm.mpfr > log.Rgemm.mpfr.2
OMP_NUM_THREADS=1 ./Rgemm.mpfr > log.Rgemm.mpfr.1

exit
