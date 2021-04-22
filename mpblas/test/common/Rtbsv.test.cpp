/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
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
#include <mpblas.h>
#include <blas.h>
#include <mplapack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N -2
#define MAX_N  7
#define MAX_K  7
#define MIN_LDA -2
#define MAX_LDA 7
#define MIN_INCX -2
#define MAX_INCX 3
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rtbsv_test2(const char *uplo, const char *trans, const char *diag)
{
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = MIN_N; n <= MAX_N; n++) {
	    for (int k = 0; k <= MAX_K; k++) {
		for (int lda = k + 1; lda <= MAX_LDA; lda++) {
#if defined VERBOSE_TEST
		    printf("#n is %d, lda is %d, incx is %d ", n, lda, incx);
		    printf("uplo is %s trans is %s, diag is %s \n", uplo, trans, diag);
#endif
		    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
		    REAL_REF *x_ref = new REAL_REF[veclen(n, incx)];
		    REAL *A = new REAL[matlen(lda, n)];
		    REAL *x = new REAL[veclen(n, incx)];

		    for (int iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(A_ref, A, matlen(lda, n));
			set_random_vector(x_ref, x, veclen(n, incx));

			mplapack_errno = 0; blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			dtbsv_f77(uplo, trans, diag, &n, &k, A_ref, &lda, x_ref, &incx);
			mplapack_errno1 = blas_errno;
#else
			Rtbsv(uplo, trans, diag, n, k, A_ref, lda, x_ref, incx);
			mplapack_errno1 = mplapack_errno;
#endif
			Rtbsv(uplo, trans, diag, n, k, A, lda, x, incx);
			mplapack_errno2 = mplapack_errno;
#if defined VERBOSE_TEST
			printf("errno: mplapack %d, ref %d\n", mplapack_errno1, mplapack_errno2);
#endif
			if (mplapack_errno1 != mplapack_errno2) {
#if defined VERBOSE_TEST
			    printf("error in Mxerbla!!\n");
#endif
			    errorflag = TRUE;
			}
			REAL_REF diff = infnorm(x_ref, x, veclen(n, incx), 1);
			if (diff > EPSILON10) {
#if defined VERBOSE_TEST
			    printf("error: "); printnum(diff); printf("\n");
#endif
			    errorflag = TRUE;
			}
			if (maxdiff < diff)
			    maxdiff = diff;
		    }
		    delete[]A_ref;
		    delete[]x_ref;
		    delete[]x;
		    delete[]A;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Rtbsv failed ***\n");
	exit(1);
    } else {
        printf("maxerror: "); printnum(maxdiff); printf("\n");
    }
}

void Rtbsv_test()
{
    Rtbsv_test2("U", "N", "U");
    Rtbsv_test2("U", "N", "N");
    Rtbsv_test2("U", "T", "U");
    Rtbsv_test2("U", "T", "N");
    Rtbsv_test2("U", "C", "U");
    Rtbsv_test2("U", "C", "N");
    Rtbsv_test2("L", "N", "U");
    Rtbsv_test2("L", "N", "N");
    Rtbsv_test2("L", "T", "U");
    Rtbsv_test2("L", "T", "N");
    Rtbsv_test2("L", "C", "U");
    Rtbsv_test2("L", "C", "N");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rtbsv start ***\n");
    Rtbsv_test();
    printf("*** Testing Rtbsv successful ***\n");
    return (0);
}
