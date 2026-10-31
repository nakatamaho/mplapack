// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_dd_cuda.h>
#include "../../../../mplapack/optimized/common/Rsolver_cpu_fallback.h"
namespace {
#define Rtrsm solver_Rtrsm_cpu
#define Rgetrs solver_reference_cpu
#include "../../../../mplapack/reference/Rgetrs.cpp"
#undef Rgetrs
#undef Rtrsm
}
void Rgetrs(const char *trans, mplapackint const n, mplapackint const nrhs, dd_real *a, mplapackint const lda, mplapackint *ipiv, dd_real *b, mplapackint const ldb, mplapackint &info) {
    if (!Rgetrs_dd_cuda(trans,n,nrhs,a,lda,ipiv,b,ldb,info)) solver_reference_cpu(trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}
