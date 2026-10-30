// Runtime complex policy, on the same resident scheduler as fixed and real.
#include "Cfactor_device_body.cuh"
namespace mplapack_mpfr_cuda {
using namespace mplapack_blas_dev;
#define INSTANTIATE(NL) template int factor_complex<factor_runtime_tag<NL>>(const factor_shape &,long,factor_runtime_complex<NL> *,const factor_runtime_complex<NL> &,long *,long &,factor_times &);
MPLAPACK_MPFR_CUDA_RT_TIERS(INSTANTIATE)
#undef INSTANTIATE
}
