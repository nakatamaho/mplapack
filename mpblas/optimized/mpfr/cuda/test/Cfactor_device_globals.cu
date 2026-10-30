// Standalone device tests link the CPU oracle, not the CUDA BLAS library.
#include "../Cblas_mpfr_dev_types.h"
namespace mplapack_blas_dev {
long mplapack_mpfr_rt_prec_h=0;
__device__ long mplapack_mpfr_rt_prec_d=0;
}
