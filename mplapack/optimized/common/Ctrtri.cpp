// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
#include <mpblas.h>
#include <mplapack.h>
#ifdef MPLAPACK_ACCELERATOR_LIBRARY
#include "Rinverse_cpu_blas.h"
#endif
#define Ctrtri Ctrtri_reference_cpu
#include "../../reference/Ctrtri.cpp"
#undef Ctrtri
#include "Rinverse_cpu_common.h"
#ifdef MPLAPACK_TRTRI_CPU_DISPATCH
static void trtri_cpu(const char *uplo,const char *diag,INTEGER n,COMPLEX *a,INTEGER lda,INTEGER &info)
#else
void Ctrtri(const char *uplo,const char *diag,INTEGER const n,COMPLEX *a,INTEGER const lda,INTEGER &info)
#endif
{
    const bool upper=Mlsame(uplo,"U");
    if((!upper && !Mlsame(uplo,"L")) || (!Mlsame(diag,"N") && !Mlsame(diag,"U")) || n<=0 || lda<std::max(INTEGER(1),n)) {
        Ctrtri_reference_cpu(uplo,diag,n,a,lda,info); return;
    }
    const INTEGER nb=iMlaenv(1,"Ctrtri",CHAR2(uplo,diag),n,-1,-1,-1);
    mplapack_inverse::tri(upper,diag,n,a,lda,nb,info);
}
