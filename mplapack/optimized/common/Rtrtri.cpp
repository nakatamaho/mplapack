// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
#include <mpblas.h>
#include <mplapack.h>
#ifdef MPLAPACK_ACCELERATOR_LIBRARY
#include "Rinverse_cpu_blas.h"
#endif
#define Rtrtri Rtrtri_reference_cpu
#include "../../reference/Rtrtri.cpp"
#undef Rtrtri
#include "Rinverse_cpu_common.h"
#ifdef MPLAPACK_TRTRI_CPU_DISPATCH
static void trtri_cpu(const char *uplo,const char *diag,INTEGER n,REAL *a,INTEGER lda,INTEGER &info)
#else
void Rtrtri(const char *uplo,const char *diag,INTEGER const n,REAL *a,INTEGER const lda,INTEGER &info)
#endif
{
    const bool upper=Mlsame(uplo,"U");
    if((!upper && !Mlsame(uplo,"L")) || (!Mlsame(diag,"N") && !Mlsame(diag,"U")) || n<=0 || lda<std::max(INTEGER(1),n)) {
        Rtrtri_reference_cpu(uplo,diag,n,a,lda,info); return;
    }
    const INTEGER nb=iMlaenv(1,"Rtrtri",CHAR2(uplo,diag),n,-1,-1,-1);
    mplapack_inverse::tri(upper,diag,n,a,lda,nb,info);
}
