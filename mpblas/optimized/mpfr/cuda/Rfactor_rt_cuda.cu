// Reuse the existing cu_mpfr arithmetic, pack and arena implementation without
// defining a second copy of its public GEMM/BLAS helpers.
#define MPLAPACK_MPFR_FACTOR_RT_ONLY
#include "Rgemm_rt_cuda.cu"
