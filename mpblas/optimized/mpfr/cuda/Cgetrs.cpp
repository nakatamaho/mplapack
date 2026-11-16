// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include "../../../../mplapack/optimized/common/Rsolver_cpu_fallback.h"
namespace {
#define Ctrsm solver_Ctrsm_cpu
#define Cgetrs solver_reference_cpu
#include "../../../../mplapack/reference/Cgetrs.cpp"
#undef Cgetrs
#undef Ctrsm
}
void Cgetrs(const char *trans, mplapackint const n, mplapackint const nrhs, mpc_class *a, mplapackint const lda, mplapackint *ipiv, mpc_class *b, mplapackint const ldb, mplapackint &info) {
    if (!Cgetrs_mpfr_cuda(trans,n,nrhs,a,lda,ipiv,b,ldb,info)) solver_reference_cpu(trans,n,nrhs,a,lda,ipiv,b,ldb,info);
}
