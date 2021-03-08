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
#define MAX_N 10
#define MIN_LDA -2
#define MAX_LDA 10
#define MIN_INCX -2
#define MAX_INCX 2
#define MIN_INCY -2
#define MAX_INCY 2
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rsymv_test3(const char *uplo, REAL_REF alpha_ref, REAL_REF beta_ref, REAL alpha, REAL beta)
{
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	    for (int n = MIN_N; n < MAX_N; n++) {
		for (int lda = max(1, n); lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
		    printf("#n is %d, lda is %d incx %d incy %d uplo %s\n", n, lda, incx, incy, uplo);
#endif
		    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
		    REAL_REF *x_ref = new REAL_REF[veclen(n, incx)];
		    REAL_REF *y_ref = new REAL_REF[veclen(n, incy)];
		    REAL *A = new REAL[matlen(lda, n)];
		    REAL *x = new REAL[veclen(n, incx)];
		    REAL *y = new REAL[veclen(n, incy)];

		    for (int iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(A_ref, A, matlen(lda, n));
			set_random_vector(x_ref, x, veclen(n, incx));
			set_random_vector(y_ref, y, veclen(n, incy));

			mplapack_errno = 0; blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			dsymv_f77(uplo, &n, &alpha_ref, A_ref, &lda, x_ref, &incx, &beta_ref, y_ref, &incy);
			mplapack_errno1 = blas_errno;
#else
			Rsymv(uplo, n, alpha_ref, A_ref, lda, x_ref, incx, beta_ref, y_ref, incy);
			mplapack_errno1 = mplapack_errno;
#endif
			Rsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy);
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
			REAL_REF diff = infnorm(y_ref, y, veclen(n, incy), 1);
			if (diff > EPSILON) {
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
		    delete[]y_ref;
		    delete[]x;
		    delete[]y;
		    delete[]A;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Rsymv failed ***\n");
	exit(1);
    }
}

void Rsymv_test2(const char *uplo)
{
    REAL_REF alpha_ref, beta_ref;
    REAL alpha, beta;

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0; beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0; beta = 0.0;
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0; beta = 1.0;
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=*;
    alpha_ref = 0.0;
    alpha = 0.0;
    set_random_number(beta_ref, beta);
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0;
    alpha = 1.0;
    set_random_number(beta_ref, beta);
    Rsymv_test3(uplo, alpha_ref, beta_ref, alpha, beta);
}

void Rsymv_test()
{
    Rsymv_test2("U");
    Rsymv_test2("L");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rsymv start ***\n");
    Rsymv_test();
    printf("*** Testing Rsymv successful ***\n");
    return (0);
}
