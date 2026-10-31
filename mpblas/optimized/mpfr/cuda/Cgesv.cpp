// CUDA dispatcher; only the entry computes on the accelerator.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#define MPLAPACK_SOLVER_CPU_COMPLEX
#include "../../../../mplapack/optimized/common/Rsolver_cpu_factor_fallback.h"
#undef MPLAPACK_SOLVER_CPU_COMPLEX
namespace {
#define Cgetrf solver_Cgetrf_cpu
#define Cgetrs solver_Cgetrs_cpu
#define Cgesv solver_reference_cpu
#include "../../../../mplapack/reference/Cgesv.cpp"
#undef Cgesv
#undef Cgetrf
#undef Cgetrs
}
void Cgesv(mplapackint const n, mplapackint const nrhs, mpc_class *a, mplapackint const lda, mplapackint *ipiv, mpc_class *b, mplapackint const ldb, mplapackint &info) {
    if (!Cgesv_mpfr_cuda(n,nrhs,a,lda,ipiv,b,ldb,info)) solver_reference_cpu(n,nrhs,a,lda,ipiv,b,ldb,info);
}
