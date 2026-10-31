// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include "../../../../mplapack/optimized/common/Rsolver_cpu_factor_fallback.h"
namespace {
#define Rgetrf solver_Rgetrf_cpu
#define Rgetrs solver_Rgetrs_cpu
#define Rgesv solver_reference_cpu
#include "../../../../mplapack/reference/Rgesv.cpp"
#undef Rgesv
#undef Rgetrf
#undef Rgetrs
}
void Rgesv(mplapackint const n, mplapackint const nrhs, mpfr_class *a, mplapackint const lda, mplapackint *ipiv, mpfr_class *b, mplapackint const ldb, mplapackint &info) {
    if (!Rgesv_mpfr_cuda(n,nrhs,a,lda,ipiv,b,ldb,info)) solver_reference_cpu(n,nrhs,a,lda,ipiv,b,ldb,info);
}
