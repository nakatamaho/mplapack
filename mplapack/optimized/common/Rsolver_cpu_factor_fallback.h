#ifndef MPLAPACK_SOLVER_CPU_FACTOR_FALLBACK_H
#define MPLAPACK_SOLVER_CPU_FACTOR_FALLBACK_H
#include "Rsolver_cpu_fallback.h"
// Bind reference factors/solves to CPU BLAS only for the new drivers' fallback.
// Existing public factor dispatchers/kernels and schedules are untouched.
#ifdef MPLAPACK_SOLVER_CPU_COMPLEX
extern decltype(iCamax) iCamax_reference_cpu;
namespace {
#define Ctrsm solver_Ctrsm_cpu
#define Cgemm solver_Cgemm_cpu
#ifdef MPLAPACK_SOLVER_CPU_CHOLESKY
#define Cherk solver_Cherk_cpu
#define Cpotrf2 solver_Cpotrf2_cpu
#define Cpotrf solver_Cpotrf_cpu
#define Cpotrs solver_Cpotrs_cpu
#include "../../reference/Cpotrf2.cpp"
#include "../../reference/Cpotrf.cpp"
#include "../../reference/Cpotrs.cpp"
#undef Cherk
#undef Cpotrf2
#undef Cpotrf
#undef Cpotrs
#else
#define Cswap solver_Cswap_cpu
#define Cscal solver_Cscal_cpu
#define iCamax iCamax_reference_cpu
#define Cgetrf2 solver_Cgetrf2_cpu
#define Cgetrf solver_Cgetrf_cpu
#define Cgetrs solver_Cgetrs_cpu
#include "../../reference/Cgetrf2.cpp"
#include "../../reference/Cgetrf.cpp"
#include "../../reference/Cgetrs.cpp"
#undef Cswap
#undef Cscal
#undef iCamax
#undef Cgetrf2
#undef Cgetrf
#undef Cgetrs
#endif
#undef Ctrsm
#undef Cgemm
}
#else
extern decltype(iRamax) iRamax_reference_cpu;
namespace {
// Legacy Rgemm_ref skips zero B entries, unlike the current reference;
// that changes signed zeros. Include the current text under a private name.
#define Rgemm solver_Rgemm_cpu
#include "../../../mpblas/reference/Rgemm.cpp"
#undef Rgemm
#define Rtrsm solver_Rtrsm_cpu
#define Rgemm solver_Rgemm_cpu
#ifdef MPLAPACK_SOLVER_CPU_CHOLESKY
#define Rsyrk solver_Rsyrk_cpu
#define Rpotrf2 solver_Rpotrf2_cpu
#define Rpotrf solver_Rpotrf_cpu
#define Rpotrs solver_Rpotrs_cpu
#include "../../reference/Rpotrf2.cpp"
#include "../../reference/Rpotrf.cpp"
#include "../../reference/Rpotrs.cpp"
#undef Rsyrk
#undef Rpotrf2
#undef Rpotrf
#undef Rpotrs
#else
#define Rswap solver_Rswap_cpu
#define Rscal solver_Rscal_cpu
#define iRamax iRamax_reference_cpu
#define Rgetrf2 solver_Rgetrf2_cpu
#define Rgetrf solver_Rgetrf_cpu
#define Rgetrs solver_Rgetrs_cpu
#include "../../reference/Rgetrf2.cpp"
#include "../../reference/Rgetrf.cpp"
#include "../../reference/Rgetrs.cpp"
#undef Rswap
#undef Rscal
#undef iRamax
#undef Rgetrf2
#undef Rgetrf
#undef Rgetrs
#endif
#undef Rtrsm
#undef Rgemm
}
#endif
#endif
