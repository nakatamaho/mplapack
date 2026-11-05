#define MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include "../Rinverse_device_cuda.cu"
#ifndef MPLAPACK_CBLAS_PROVIDES_RT_GLOBALS
namespace mplapack_blas_dev { long mplapack_mpfr_rt_prec_h=0; }
#endif
