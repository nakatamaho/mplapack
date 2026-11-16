#include <mpblas_dd_cuda.h>
#define MPLAPACK_POTRF_CPU_DISPATCH
#include "../../../../mplapack/optimized/common/Cpotrf.cpp"
void Cpotrf(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER &info) {
    if (!Cpotrf_dd_cuda(uplo,n,a,lda,info)) potrf_cpu(uplo,n,a,lda,info);
}
