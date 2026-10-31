// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_td_cuda.h>
#include "../../../../mplapack/optimized/common/Rsolver_cpu_fallback.h"
namespace {
#define Ctrsm solver_Ctrsm_cpu
#define Cpotrs solver_reference_cpu
#include "../../../../mplapack/reference/Cpotrs.cpp"
#undef Cpotrs
#undef Ctrsm
}
void Cpotrs(const char *uplo, mplapackint const n, mplapackint const nrhs, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, mplapackint &info) {
    if (!Cpotrs_td_cuda(uplo,n,nrhs,a,lda,b,ldb,info)) solver_reference_cpu(uplo,n,nrhs,a,lda,b,ldb,info);
}
