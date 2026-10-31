#ifndef MPLAPACK_SOLVER_LEGACY_BENCH_H
#define MPLAPACK_SOLVER_LEGACY_BENCH_H
// Historical CPU outer solvers with per-call GPU BLAS. The new dispatcher
// deliberately uses CPU TRSM on rejection, so it is no longer this baseline.
// The runner also sets <P>_LAPACK_CUDA=0 to disable resident factorizations.
namespace {
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
}
#endif
