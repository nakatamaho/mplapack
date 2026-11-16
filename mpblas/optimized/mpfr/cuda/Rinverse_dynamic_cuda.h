// Internal, dynamically sized real runtime inverse engine.
#ifndef MPLAPACK_MPFR_INVERSE_DYNAMIC_CUDA_H
#define MPLAPACK_MPFR_INVERSE_DYNAMIC_CUDA_H
#include "Rinverse_kernel_cuda.h"
namespace mplapack_mpfr_inverse {
__attribute__((visibility("hidden"))) bool dynamic_real_memory_fits(const shape &,long);
__attribute__((visibility("hidden"))) int run_dynamic_real(const shape &,long,mplapack_mpfr_cuda::rt_pack &,const long *,long &,mplapack_mpfr_cuda::factor_times &);
}
#endif
