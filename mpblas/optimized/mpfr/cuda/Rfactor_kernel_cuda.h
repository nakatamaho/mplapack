// Internal MPFR device-resident real factorization interface. No MPFR C++ types.
#ifndef MPLAPACK_MPFR_FACTOR_KERNEL_CUDA_H
#define MPLAPACK_MPFR_FACTOR_KERNEL_CUDA_H
#include "Rgemm_kernel_cuda.h"
#include "Rgemm_rt_cuda.h"
#include "../../common/Rpotrf_output_common.h"
namespace mplapack_mpfr_cuda {
struct factor_shape {
    bool chol, upper;
    long m, n, nb;
    // Optional solver continuation. Old factor-only calls leave these defaults.
    long nrhs=0;
    bool factorize=true;
    char trans='N';
};
static inline bool factor_output(const factor_shape &s,long info,long i,long j) {
    return j>=s.n || !s.chol || mplapack_lapack_gpu::potrf_output(s.upper,s.n,s.nb,info,i,j);
}
struct factor_times { double upload=0, factor=0, download=0; };
__attribute__((visibility("hidden"))) size_t factor_stack_bytes(long prec,bool runtime);
__attribute__((visibility("hidden"))) bool factor_memory_fits(long m,long n,long prec,bool fixed);
template<int PB> __attribute__((visibility("hidden"))) int factor_fixed(const factor_shape &, real_t<PB> *, const real_t<PB> &, long *, long &, factor_times &);
__attribute__((visibility("hidden"))) int factor_runtime(const factor_shape &, long, rt_pack &, const rt_pack &, long *, long &, factor_times &);
}
#endif
