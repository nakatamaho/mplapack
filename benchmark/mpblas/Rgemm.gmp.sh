#!/bin/bash

OMP_NUM_THREADS=64 ./Rgemm.gmp > log.Rgemm.gmp.64
OMP_NUM_THREADS=60 ./Rgemm.gmp > log.Rgemm.gmp.60
OMP_NUM_THREADS=48 ./Rgemm.gmp > log.Rgemm.gmp.48
OMP_NUM_THREADS=32 ./Rgemm.gmp > log.Rgemm.gmp.32
OMP_NUM_THREADS=30 ./Rgemm.gmp > log.Rgemm.gmp.30
OMP_NUM_THREADS=24 ./Rgemm.gmp > log.Rgemm.gmp.24
OMP_NUM_THREADS=16 ./Rgemm.gmp > log.Rgemm.gmp.16
OMP_NUM_THREADS=8 ./Rgemm.gmp > log.Rgemm.gmp.8
OMP_NUM_THREADS=4 ./Rgemm.gmp > log.Rgemm.gmp.4
OMP_NUM_THREADS=2 ./Rgemm.gmp > log.Rgemm.gmp.2
OMP_NUM_THREADS=1 ./Rgemm.gmp > log.Rgemm.gmp.1

exit
