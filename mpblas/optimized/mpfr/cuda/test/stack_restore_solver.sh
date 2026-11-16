#!/bin/sh
exec env MPLAPACK_MPFR_CUDA_STACK_RESTORE=1 MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0 ./solver_gpu_device_test_mpfr 256 --quick
