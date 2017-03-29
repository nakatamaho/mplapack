#!/bin/bash

OMP_NUM_THREADS=64 ./Rgemm.double > log.Rgemm.double.64
OMP_NUM_THREADS=60 ./Rgemm.double > log.Rgemm.double.60
OMP_NUM_THREADS=48 ./Rgemm.double > log.Rgemm.double.48
OMP_NUM_THREADS=32 ./Rgemm.double > log.Rgemm.double.32
OMP_NUM_THREADS=30 ./Rgemm.double > log.Rgemm.double.30
OMP_NUM_THREADS=24 ./Rgemm.double > log.Rgemm.double.24
OMP_NUM_THREADS=16 ./Rgemm.double > log.Rgemm.double.16
OMP_NUM_THREADS=8 ./Rgemm.double > log.Rgemm.double.8
OMP_NUM_THREADS=4 ./Rgemm.double > log.Rgemm.double.4
OMP_NUM_THREADS=2 ./Rgemm.double > log.Rgemm.double.2
OMP_NUM_THREADS=1 ./Rgemm.double > log.Rgemm.double.1

exit
