#ifndef MPLAPACK_SOLVER_CPU_FACTOR_FALLBACK_H
#define MPLAPACK_SOLVER_CPU_FACTOR_FALLBACK_H
#include "Rsolver_cpu_fallback.h"
// Driver fallback = public factor dispatcher + CPU-only solve. In particular,
// do not replace the optimized factor/panel/GEMM by private serial copies.
// Drivers follow the public MPFR factor's signed-zero conventions.
namespace {
#ifdef MPLAPACK_SOLVER_CPU_COMPLEX
#define Ctrsm solver_Ctrsm_cpu
#ifdef MPLAPACK_SOLVER_CPU_CHOLESKY
#define Cpotrs solver_Cpotrs_cpu
#include "../../reference/Cpotrs.cpp"
#undef Cpotrs
#else
#define Cgetrs solver_Cgetrs_cpu
#include "../../reference/Cgetrs.cpp"
#undef Cgetrs
#endif
#undef Ctrsm
#else
#define Rtrsm solver_Rtrsm_cpu
#ifdef MPLAPACK_SOLVER_CPU_CHOLESKY
#define Rpotrs solver_Rpotrs_cpu
#include "../../reference/Rpotrs.cpp"
#undef Rpotrs
#else
#define Rgetrs solver_Rgetrs_cpu
#include "../../reference/Rgetrs.cpp"
#undef Rgetrs
#endif
#undef Rtrsm
#endif
}
#endif
