#ifndef MPLAPACK_SOLVER_LEGACY_BENCH_H
#define MPLAPACK_SOLVER_LEGACY_BENCH_H
// Historical CPU outer solvers with per-call GPU BLAS. The new dispatcher
// deliberately uses CPU TRSM on rejection, so it is no longer this baseline.
// Public factor CPU fallbacks are now CPU-only too. Keep private generated
// factors here so --legacy remains a per-call GPU-BLAS comparator, not a
// hybrid of the new CPU factors and old GPU solves. Benchmark code only.
namespace {
#define Rgetrf2 solver_legacy_Rgetrf2
#define Cgetrf2 solver_legacy_Cgetrf2
#define Rpotrf2 solver_legacy_Rpotrf2
#define Cpotrf2 solver_legacy_Cpotrf2
#define Rgetrf solver_legacy_Rgetrf
#define Cgetrf solver_legacy_Cgetrf
#define Rpotrf solver_legacy_Rpotrf
#define Cpotrf solver_legacy_Cpotrf
#include "../mplapack/reference/Rgetrf2.cpp"
#include "../mplapack/reference/Cgetrf2.cpp"
#include "../mplapack/reference/Rpotrf2.cpp"
#include "../mplapack/reference/Cpotrf2.cpp"
#include "../mplapack/reference/Rgetrf.cpp"
#include "../mplapack/reference/Cgetrf.cpp"
#include "../mplapack/reference/Rpotrf.cpp"
#include "../mplapack/reference/Cpotrf.cpp"
#define Rgetrs solver_legacy_Rgetrs
#define Cgetrs solver_legacy_Cgetrs
#define Rpotrs solver_legacy_Rpotrs
#define Cpotrs solver_legacy_Cpotrs
#define Rgesv solver_legacy_Rgesv
#define Cgesv solver_legacy_Cgesv
#define Rposv solver_legacy_Rposv
#define Cposv solver_legacy_Cposv
#include "../mplapack/reference/Rgetrs.cpp"
#include "../mplapack/reference/Cgetrs.cpp"
#include "../mplapack/reference/Rpotrs.cpp"
#include "../mplapack/reference/Cpotrs.cpp"
#include "../mplapack/reference/Rgesv.cpp"
#include "../mplapack/reference/Cgesv.cpp"
#include "../mplapack/reference/Rposv.cpp"
#include "../mplapack/reference/Cposv.cpp"
#undef Rgetrs
#undef Cgetrs
#undef Rpotrs
#undef Cpotrs
#undef Rgesv
#undef Cgesv
#undef Rposv
#undef Cposv
#undef Rgetrf2
#undef Cgetrf2
#undef Rpotrf2
#undef Cpotrf2
#undef Rgetrf
#undef Cgetrf
#undef Rpotrf
#undef Cpotrf
}
#endif
