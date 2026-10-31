// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_dd_cuda.h>
#define MPLAPACK_SOLVER_CPU_CHOLESKY
#include "../../../../mplapack/optimized/common/Rsolver_cpu_factor_fallback.h"
#undef MPLAPACK_SOLVER_CPU_CHOLESKY
namespace {
#define Rpotrf solver_Rpotrf_cpu
#define Rpotrs solver_Rpotrs_cpu
#define Rposv solver_reference_cpu
#include "../../../../mplapack/reference/Rposv.cpp"
#undef Rposv
#undef Rpotrf
#undef Rpotrs
}
void Rposv(const char *uplo, mplapackint const n, mplapackint const nrhs, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, mplapackint &info) {
    if (!Rposv_dd_cuda(uplo,n,nrhs,a,lda,b,ldb,info)) solver_reference_cpu(uplo,n,nrhs,a,lda,b,ldb,info);
}
