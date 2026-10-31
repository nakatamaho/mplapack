// CPU-library bindings for accelerator inverse fallbacks, including leaves.
#ifndef MPLAPACK_INVERSE_CPU_BLAS_H
#define MPLAPACK_INVERSE_CPU_BLAS_H
#include "Rfactor_cpu_blas.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_MPFR)
#define MPLAPACK_INVERSE_CPU_BLAS(R) \
extern mplapack_solver_cpu_detail::boolean_result<decltype(R)>::type R##_omp; \
extern decltype(R) R##_reference_cpu; \
namespace { template<class... A> void inverse_##R##_cpu(A... a) { \
    if(!R##_omp(a...)) R##_reference_cpu(a...); \
} }
MPLAPACK_INVERSE_CPU_BLAS(Rtrmm)
MPLAPACK_INVERSE_CPU_BLAS(Ctrmm)
MPLAPACK_INVERSE_CPU_BLAS(Rtrmv)
MPLAPACK_INVERSE_CPU_BLAS(Ctrmv)
MPLAPACK_INVERSE_CPU_BLAS(Rgemv)
MPLAPACK_INVERSE_CPU_BLAS(Cgemv)
MPLAPACK_INVERSE_CPU_BLAS(CRscal)
#undef MPLAPACK_INVERSE_CPU_BLAS
#define Rtrmm inverse_Rtrmm_cpu
#define Ctrmm inverse_Ctrmm_cpu
#define Rtrmv inverse_Rtrmv_cpu
#define Ctrmv inverse_Ctrmv_cpu
#define Rgemv inverse_Rgemv_cpu
#define Cgemv inverse_Cgemv_cpu
#define CRscal inverse_CRscal_cpu
// Rdot/Cdotc already resolve to CPU implementations in these CUDA flavors.
#define Rtrti2 inverse_Rtrti2_cpu
#define Ctrti2 inverse_Ctrti2_cpu
#define Rlauu2 inverse_Rlauu2_cpu
#define Clauu2 inverse_Clauu2_cpu
namespace {
#include "../../reference/Rtrti2.cpp"
#include "../../reference/Ctrti2.cpp"
#include "../../reference/Rlauu2.cpp"
#include "../../reference/Clauu2.cpp"
}
#endif
#endif
