#include "Rinverse_device_body.cuh"
namespace mplapack_mpfr_inverse {
#define INSTANTIATE(NL) \
template int run<factor_runtime_tag<NL>,false>(const shape &,long,value<factor_runtime_tag<NL>,false> *,value<factor_runtime_tag<NL>,false> *,const long *,long &,factor_times &); \
template int run<factor_runtime_tag<NL>,true>(const shape &,long,value<factor_runtime_tag<NL>,true> *,value<factor_runtime_tag<NL>,true> *,const long *,long &,factor_times &);
MPLAPACK_MPFR_CUDA_RT_TIERS(INSTANTIATE)
#undef INSTANTIATE
}
