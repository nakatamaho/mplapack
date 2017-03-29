/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rpotri.debug.cpp,v 1.2 2010/08/19 01:17:55 nakatamaho Exp $
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
#include <mblas.h>
#include <mlapack.h>
#include <mpack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N      1
#define MAX_N     12
#define MAX_LDA   12
#define MAX_ITER   5

REAL_REF maxdiff = 0.0;

void Rpotri_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    REAL_REF diff;
    INTEGER info;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
	    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
	    REAL *A = new REAL[matlen(lda, n)];

#if defined VERBOSE_TEST
	    printf("#n:%d lda %d\n", n, lda);
#endif
	    j = 0;
	    while (j < MAX_ITER) {
		set_random_symmmat_cond(A_ref, A, lda, n, 2);
//		set_random_psdmat(A_ref, A, lda, n);
//numerical error measure: first do inversion
#if defined ___MPACK_BUILD_WITH_MPFR___
		dpotf2_f77(uplo, &n, A_ref, &lda, &info_ref);
		dpotri_f77(uplo, &n, A_ref, &lda, &info_ref);
#else
		Rpotrf(uplo, n, A_ref, lda, &info_ref);
		Rpotri(uplo, n, A_ref, lda, &info_ref);
#endif
	      for (int p = 0; p < matlen(lda, n); p++ ) {
#if defined ___MPACK_BUILD_WITH_MPFR___
		 A[p] = A_ref[p];
#elif defined ___MPACK_BUILD_WITH_GMP___
		 A[p] = cast2mpf_class(A_ref[p]);
#elif defined ___MPACK_BUILD_WITH_QD___
		 A[p] = cast2qd_real(A_ref[p]);
#elif defined ___MPACK_BUILD_WITH_DD___
		 A[p] = cast2dd_real(A_ref[p]);
#elif defined ___MPACK_BUILD_WITH_DOUBLE___
		 A[p] = cast2double(A_ref[p]);
#elif defined ___MPACK_BUILD_WITH___FLOAT128___
		 A[p] = cast2__float128(A_ref[p]);
#endif
	      }
//doing inversion twice.
#if defined ___MPACK_BUILD_WITH_MPFR___
		dpotf2_f77(uplo, &n, A_ref, &lda, &info_ref);
		dpotri_f77(uplo, &n, A_ref, &lda, &info_ref);
#else
		Rpotrf(uplo, n, A_ref, lda, &info_ref);
		Rpotri(uplo, n, A_ref, lda, &info_ref);
#endif
		Rpotrf(uplo, n, A, lda, &info);
		Rpotri(uplo, n, A, lda, &info);

		if (info < 0) {
#if defined VERBOSE_TEST
		    printf("info %d error\n", -(int) info);
#endif
		}
		if (info_ref != info) {
		    printf("info differ! %d, %d\n", (int) info_ref, (int) info);
		    errorflag = TRUE;
		}
		diff = infnorm_mat(A_ref, A, n, n, lda);
		if (diff > EPSILON6) {
		    printf("error: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
                printf("A_ref = "); printmat (n,n,A_ref,lda); printf("\n");
                printf("A ="); printmat (n,n,A,lda); printf("\n");
		}
	        if (maxdiff < diff)  maxdiff = diff;
#if defined VERBOSE_TEST
	        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		j++;
	    }
	    delete[]A_ref;
	    delete[]A;
	}
	if (errorflag == TRUE) {
            printf("*** Testing Rpotri failed ***\n");
	    exit(1);
	}
    }
}
void Rpotri_test()
{
     Rpotri_test2("U");
     Rpotri_test2("L");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rpotri start ***\n");
    Rpotri_test();
    printf("*** Testing Rpotri successful ***\n");
    return (0);
}
