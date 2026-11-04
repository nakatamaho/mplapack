#!/bin/sh
set -eu
for p in 512 200; do
    MPFRXX_DEFAULT_PRECISION_BITS=$p OMP_NUM_THREADS=1 MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0 ./inverse_gpu_host_test_mpfr
done
MPFRXX_DEFAULT_PRECISION_BITS=512 OMP_NUM_THREADS=1 MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1 ./inverse_gpu_host_test_mpfr
MPFRXX_DEFAULT_PRECISION_BITS=8192 OMP_NUM_THREADS=1 ./inverse_gpu_host_test_mpfr --real-blocked
MPFRXX_DEFAULT_PRECISION_BITS=16384 OMP_NUM_THREADS=1 ./inverse_gpu_host_test_mpfr --real-small
