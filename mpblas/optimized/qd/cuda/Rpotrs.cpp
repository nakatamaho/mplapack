// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_qd_cuda.h>
#include "../../../../mplapack/optimized/common/Rsolver_cpu_fallback.h"
namespace {
#define Rtrsm solver_Rtrsm_cpu
#define Rpotrs solver_reference_cpu
#include "../../../../mplapack/reference/Rpotrs.cpp"
#undef Rpotrs
#undef Rtrsm
}
void Rpotrs(const char *uplo, mplapackint const n, mplapackint const nrhs, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb, mplapackint &info) {
    if (!Rpotrs_qd_cuda(uplo,n,nrhs,a,lda,b,ldb,info)) solver_reference_cpu(uplo,n,nrhs,a,lda,b,ldb,info);
}
