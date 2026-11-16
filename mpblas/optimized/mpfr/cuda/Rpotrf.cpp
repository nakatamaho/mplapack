#include <mpblas_mpfr_cuda.h>
#define MPLAPACK_POTRF_CPU_DISPATCH
#include "../../../../mplapack/optimized/common/Rpotrf.cpp"
void Rpotrf(const char *uplo,INTEGER const n,REAL *a,INTEGER const lda,INTEGER &info) {
    if(!Rpotrf_mpfr_cuda(uplo,n,a,lda,info)) potrf_cpu(uplo,n,a,lda,info);
}
