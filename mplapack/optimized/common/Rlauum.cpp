// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
#include <mpblas.h>
#include <mplapack.h>
#ifdef MPLAPACK_ACCELERATOR_LIBRARY
#include "Rinverse_cpu_blas.h"
#endif
#define Rlauum Rlauum_reference_cpu
#include "../../reference/Rlauum.cpp"
#undef Rlauum
#include "Rinverse_cpu_common.h"
#ifdef MPLAPACK_LAUUM_CPU_DISPATCH
static void lauum_cpu(const char *uplo,INTEGER n,REAL *a,INTEGER lda,INTEGER &info)
#else
void Rlauum(const char *uplo,INTEGER const n,REAL *a,INTEGER const lda,INTEGER &info)
#endif
{
    const bool upper=Mlsame(uplo,"U");
    if((!upper && !Mlsame(uplo,"L")) || n<=0 || lda<std::max(INTEGER(1),n)) {
        Rlauum_reference_cpu(uplo,n,a,lda,info); return;
    }
    const INTEGER nb=iMlaenv(1,"Rlauum",uplo,n,-1,-1,-1);
    mplapack_inverse::lauum(upper,n,a,lda,nb,info);
}
