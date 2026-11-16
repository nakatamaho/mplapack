// Standalone MPFR512 resident benchmark: CPU oracle library + production jobs.
// No public accelerator entry is invoked or timed. BSD-2-Clause.
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#ifndef MPLAPACK_BUILD_WITH_MPFR
#define MPLAPACK_BUILD_WITH_MPFR
#endif
// nvcc's glibc feature workaround hides the GNU clock-lock declarations.
// This standalone TU uses only ordinary std::mutex, not timed locks.
#include <cstddef>
#ifdef __GLIBCXX__
#undef _GLIBCXX_USE_PTHREAD_COND_CLOCKWAIT
#undef _GLIBCXX_USE_PTHREAD_MUTEX_CLOCKLOCK
#endif
#define KERNEL_BENCH_MPFR_RESIDENT
#include "Rkernel_bench.cpp"
namespace mplapack_blas_dev {
long mplapack_mpfr_rt_prec_h=512;
__device__ long mplapack_mpfr_rt_prec_d=512;
}
namespace mplapack_mpfr_inverse {
thread_local bool failed_host=false;
__device__ bool failed_device[8192];
}
