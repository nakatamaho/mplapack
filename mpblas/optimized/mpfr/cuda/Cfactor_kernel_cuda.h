// Internal MPFR complex factorization interface; no host MPFR C++ types.
#ifndef MPLAPACK_MPFR_CFACTOR_KERNEL_CUDA_H
#define MPLAPACK_MPFR_CFACTOR_KERNEL_CUDA_H
#include "Rfactor_kernel_cuda.h"
#include "Cfactor_runtime_types.h"
namespace mplapack_mpfr_cuda {
template<class Tag> __attribute__((visibility("hidden"))) int factor_complex(
    const factor_shape &,long,typename mplapack_blas_dev::cgpu_types<Tag>::DC *,
    const typename mplapack_blas_dev::cgpu_types<Tag>::DC &,long *,long &,factor_times &);
__attribute__((visibility("hidden"))) bool complex_factor_memory_fits(long,long,long,bool);
}
#endif
