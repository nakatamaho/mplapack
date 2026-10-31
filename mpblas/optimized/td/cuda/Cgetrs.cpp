// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_td_cuda.h>
#include "../../../../mplapack/optimized/common/Rsolver_cpu_fallback.h"
namespace {
#define Ctrsm solver_Ctrsm_cpu
#define Cgetrs solver_reference_cpu
#include "../../../../mplapack/reference/Cgetrs.cpp"
#undef Cgetrs
#undef Ctrsm
}
void Cgetrs(const char *trans, mplapackint const n, mplapackint const nrhs, td_complex *a, mplapackint const lda, mplapackint *ipiv, td_complex *b, mplapackint const ldb, mplapackint &info) {
    if (!Cgetrs_td_cuda(trans,n,nrhs,a,lda,ipiv,b,ldb,info)) solver_reference_cpu(trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}
