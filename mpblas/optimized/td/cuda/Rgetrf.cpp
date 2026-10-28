#include <mpblas_td_cuda.h>
#define MPLAPACK_GETRF_CPU_DISPATCH
#include "../../../../mplapack/optimized/common/Rgetrf.cpp"
void Rgetrf(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info) {
    if (!Rgetrf_td_cuda(m,n,a,lda,ipiv,info)) getrf_cpu(m,n,a,lda,ipiv,info);
}
