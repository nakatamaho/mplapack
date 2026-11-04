// Rkernel_bench.td_cuda_kernel: the td GPU kernels on device-resident
// matrices (see Rkernel_bench.cpp).  Compile with -fmad=false and
// -Xcompiler=-ffp-contract=off, as the library's .cu files.
#ifndef MPLAPACK_BUILD_WITH_TD
#define MPLAPACK_BUILD_WITH_TD
#endif
#define KERNEL_BENCH_BF_CUDA
#include "Rkernel_bench.cpp"
