// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_td_cuda.h>
#define MPLAPACK_SOLVER_CPU_COMPLEX
#define MPLAPACK_SOLVER_CPU_CHOLESKY
#include "../../../../mplapack/optimized/common/Rsolver_cpu_factor_fallback.h"
#undef MPLAPACK_SOLVER_CPU_COMPLEX
#undef MPLAPACK_SOLVER_CPU_CHOLESKY
namespace {
#define Cpotrf solver_Cpotrf_cpu
#define Cpotrs solver_Cpotrs_cpu
#define Cposv solver_reference_cpu
#include "../../../../mplapack/reference/Cposv.cpp"
#undef Cposv
#undef Cpotrf
#undef Cpotrs
}
void Cposv(const char *uplo, mplapackint const n, mplapackint const nrhs, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, mplapackint &info) {
    if (!Cposv_td_cuda(uplo,n,nrhs,a,lda,b,ldb,info)) solver_reference_cpu(uplo,n,nrhs,a,lda,b,ldb,info);
}
