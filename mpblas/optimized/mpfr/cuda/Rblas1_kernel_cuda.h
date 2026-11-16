#ifndef MPLAPACK_MPFR_BLAS1_KERNEL_CUDA_H
#define MPLAPACK_MPFR_BLAS1_KERNEL_CUDA_H
#include "Rinverse_types.h"
#include "Rinverse_dynamic_cuda.h"
#define MPLAPACK_BLAS1_HD MPLAPACK_DEV_HD
#include "../../common/Rblas1_reduction_common.h"
namespace mplapack_mpfr_inverse {
// Workspace states of the one-pass iamax: 2 * factor_pool (8192 = 64 threads *
// 128 blocks, Rfactor_device_body.cuh), scan partials plus combine partials.
constexpr long blas1_iamax_states=2*8192;
__attribute__((visibility("hidden"))) bool dynamic_blas1_memory_fits(long,long);
__attribute__((visibility("hidden"))) int run_dynamic_blas1(mplapack_blas1::operation,long,long,mplapack_mpfr_cuda::rt_pack &,long &);
template<class Tag,bool Complex> int run_blas1(mplapack_blas1::operation,long,long,
    value<Tag,Complex> *,value<Tag,Complex> *,value<Tag,false>,mplapack_blas1::state<value<Tag,false>> &);
}
#endif
