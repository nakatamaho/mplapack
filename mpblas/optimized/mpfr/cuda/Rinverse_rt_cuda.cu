#include "Rinverse_device_body.cuh"
#include "Rblas1_device_body.cuh"
namespace mplapack_mpfr_inverse {
#define INSTANTIATE(NL) \
template int run<factor_runtime_tag<NL>,false>(const shape &,long,value<factor_runtime_tag<NL>,false> *,value<factor_runtime_tag<NL>,false> *,const long *,long &,factor_times &); \
template int run<factor_runtime_tag<NL>,true>(const shape &,long,value<factor_runtime_tag<NL>,true> *,value<factor_runtime_tag<NL>,true> *,const long *,long &,factor_times &);
MPLAPACK_MPFR_CUDA_RT_TIERS(INSTANTIATE)
#undef INSTANTIATE
#define INSTANTIATE_BLAS1(NL) \
template int run_blas1<factor_runtime_tag<NL>,false>(mplapack_blas1::operation,long,long,value<factor_runtime_tag<NL>,false> *,value<factor_runtime_tag<NL>,false> *,value<factor_runtime_tag<NL>,false>,mplapack_blas1::state<value<factor_runtime_tag<NL>,false>> &); \
template int run_blas1<factor_runtime_tag<NL>,true>(mplapack_blas1::operation,long,long,value<factor_runtime_tag<NL>,true> *,value<factor_runtime_tag<NL>,true> *,value<factor_runtime_tag<NL>,false>,mplapack_blas1::state<value<factor_runtime_tag<NL>,false>> &);
MPLAPACK_MPFR_CUDA_RT_TIERS(INSTANTIATE_BLAS1)
#undef INSTANTIATE_BLAS1
}
