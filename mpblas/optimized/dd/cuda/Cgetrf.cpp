#include <mpblas_dd_cuda.h>
#define MPLAPACK_GETRF_CPU_DISPATCH
#include "../../../../mplapack/optimized/common/Cgetrf.cpp"
void Cgetrf(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info) {
    if (!Cgetrf_dd_cuda(m,n,a,lda,ipiv,info)) getrf_cpu(m,n,a,lda,ipiv,info);
}
