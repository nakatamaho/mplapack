#!/bin/bash

OMP_NUM_THREADS=64 ./Rgemm.qd > log.Rgemm.qd.64
OMP_NUM_THREADS=60 ./Rgemm.qd > log.Rgemm.qd.60
OMP_NUM_THREADS=48 ./Rgemm.qd > log.Rgemm.qd.48
OMP_NUM_THREADS=32 ./Rgemm.qd > log.Rgemm.qd.32
OMP_NUM_THREADS=30 ./Rgemm.qd > log.Rgemm.qd.30
OMP_NUM_THREADS=24 ./Rgemm.qd > log.Rgemm.qd.24
OMP_NUM_THREADS=16 ./Rgemm.qd > log.Rgemm.qd.16
OMP_NUM_THREADS=8 ./Rgemm.qd > log.Rgemm.qd.8
OMP_NUM_THREADS=4 ./Rgemm.qd > log.Rgemm.qd.4
OMP_NUM_THREADS=2 ./Rgemm.qd > log.Rgemm.qd.2
OMP_NUM_THREADS=1 ./Rgemm.qd > log.Rgemm.qd.1

exit
