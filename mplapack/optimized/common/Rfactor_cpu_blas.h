// CPU BLAS bindings for the CPU side of CUDA factor dispatchers only.
// Reuse the CPU library's GEMM dispatcher, not the single-threaded reference
// copy or the legacy zero-skipping Rgemm_ref. No process-wide environment
// changes: unrelated/concurrent BLAS calls retain their CUDA policy.
#ifndef MPLAPACK_FACTOR_CPU_BLAS_H
#define MPLAPACK_FACTOR_CPU_BLAS_H
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_MPFR)
#include <mpblas.h>
#include <mplapack.h>
#include "Rsolver_cpu_fallback.h"
extern decltype(iRamax) iRamax_reference_cpu;
extern decltype(iCamax) iCamax_reference_cpu;
#if defined(MPLAPACK_BUILD_WITH_MPFR)
// As in the MPFR bridges, inline the private dependency graph rather than
// exporting new weak implementation instantiations from the shared library.
static decltype(Rgemm) factor_Rgemm_cpu __attribute__((flatten));
#else
static decltype(Rgemm) factor_Rgemm_cpu;
#endif
#define Rgemm factor_Rgemm_cpu
#if defined(MPLAPACK_BUILD_WITH_DD)
#include "../../../mpblas/optimized/dd/Rgemm.cpp"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include "../../../mpblas/optimized/td/Rgemm.cpp"
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include "../../../mpblas/optimized/qd/Rgemm.cpp"
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
// Autotools' MPFR CUDA library lacks Rgemm_blocked_omp (CMake has it).
// Keep this dependency private, and bind optional Winograd leaves to our
// CPU GEMM as well, rather than accidentally dispatching them to CUDA.
__attribute__((flatten)) static bool factor_Rgemm_blocked_omp(bool, bool, INTEGER, INTEGER, INTEGER,
    REAL, REAL *, INTEGER, REAL *, INTEGER, REAL, REAL *, INTEGER);
#define Rgemm_blocked_omp factor_Rgemm_blocked_omp
#include "../../../mpblas/optimized/common/blas/Rgemm_blocked_omp.cpp"
#include "../../../mpblas/optimized/mpfr/Rgemm.cpp"
// This newly needed inline constructor is not part of the shared ABI.
template __attribute__((visibility("hidden"))) mpfrxx::mpfr_class::mpfr_class(int);
#undef Rgemm_blocked_omp
#endif
#undef Rgemm
#define Rgemm factor_Rgemm_cpu
#define Cgemm solver_Cgemm_cpu
#define Rtrsm solver_Rtrsm_cpu
#define Ctrsm solver_Ctrsm_cpu
#define Rsyrk solver_Rsyrk_cpu
#define Cherk solver_Cherk_cpu
#define Rscal solver_Rscal_cpu
#define Cscal solver_Cscal_cpu
#define Rswap solver_Rswap_cpu
#define Cswap solver_Cswap_cpu
#define iRamax iRamax_reference_cpu
#define iCamax iCamax_reference_cpu
// CPU Cholesky uses the generated recursive diagonal panel as before, but
// every BLAS call of that recursion must have the same CPU binding.
namespace {
#define Rpotrf2 factor_Rpotrf2_cpu
#define Cpotrf2 factor_Cpotrf2_cpu
#include "../../reference/Rpotrf2.cpp"
#include "../../reference/Cpotrf2.cpp"
#undef Rpotrf2
#undef Cpotrf2
}
#define Rpotrf2 factor_Rpotrf2_cpu
#define Cpotrf2 factor_Cpotrf2_cpu
#endif // CUDA factor backends
#endif
