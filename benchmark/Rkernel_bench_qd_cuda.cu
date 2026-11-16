// Rkernel_bench.qd_cuda_kernel: the qd GPU kernels on device-resident
// matrices (see Rkernel_bench.cpp).  Compile with -fmad=false and
// -Xcompiler=-ffp-contract=off, as the library's .cu files.
#ifndef MPLAPACK_BUILD_WITH_QD
#define MPLAPACK_BUILD_WITH_QD
#endif
#define KERNEL_BENCH_BF_CUDA
#include "Rkernel_bench.cpp"
