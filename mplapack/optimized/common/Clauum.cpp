// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
#include <mpblas.h>
#include <mplapack.h>
#ifdef MPLAPACK_ACCELERATOR_LIBRARY
#include "Rinverse_cpu_blas.h"
#endif
#define Clauum Clauum_reference_cpu
#include "../../reference/Clauum.cpp"
#undef Clauum
#include "Rinverse_cpu_common.h"
#ifdef MPLAPACK_LAUUM_CPU_DISPATCH
static void lauum_cpu(const char *uplo,INTEGER n,COMPLEX *a,INTEGER lda,INTEGER &info)
#else
void Clauum(const char *uplo,INTEGER const n,COMPLEX *a,INTEGER const lda,INTEGER &info)
#endif
{
    const bool upper=Mlsame(uplo,"U");
    if((!upper && !Mlsame(uplo,"L")) || n<=0 || lda<std::max(INTEGER(1),n)) {
        Clauum_reference_cpu(uplo,n,a,lda,info); return;
    }
    const INTEGER nb=iMlaenv(1,"Clauum",uplo,n,-1,-1,-1);
    mplapack_inverse::lauum(upper,n,a,lda,nb,info);
}
