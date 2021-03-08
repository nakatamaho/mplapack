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

#define MIN_N -1
#define MAX_N  8
#define MAX_K  8
#define MIN_LDA 1
#define MAX_LDA 8
#define MIN_INCX -2
#define MAX_INCX 2
#define MIN_INCY -2
#define MAX_INCY 2
#define MAX_ITER 2

REAL_REF maxdiff = 0.0;

void Rsbmv_test3(const char *uplo, REAL_REF alpha_ref, REAL_REF beta_ref, REAL alpha, REAL beta)
{
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	    for (int n = MIN_N; n < MAX_N; n++) {
		for (int k = 0; k < MAX_K; k++) {
		    for (int lda = k + 1; lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
			printf("#n is %d, lda is %d, incx is %d, incy is %d ", n, lda, incx, incy);
			printf("uplo is %s \n", uplo);
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
			    mplapack_errno = 0;  blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			    dsbmv_f77(uplo, &n, &k, &alpha_ref, A_ref, &lda, x_ref, &incx, &beta_ref, y_ref, &incy);
			    mplapack_errno1 = blas_errno;
#else
			    Rsbmv(uplo, n, k, alpha_ref, A_ref, lda, x_ref, incx, beta_ref, y_ref, incy);
			    mplapack_errno1 = mplapack_errno;
#endif
			    Rsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
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
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Rsbmv failed ***\n");
	exit(1);
    }
}

void Rsbmv_test2(const char *uplo)
{
    REAL_REF alpha_ref, beta_ref;
    REAL alpha, beta;

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0; beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0; beta = 0.0;
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0; beta = 1.0;
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=*;
    alpha_ref = 0.0; alpha = 0.0;
    set_random_number(beta_ref, beta);
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0; alpha = 1.0;
    set_random_number(beta_ref, beta);
    Rsbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);
}

void Rsbmv_test()
{
    Rsbmv_test2("U");
    Rsbmv_test2("L");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rsbmv start ***\n");
    Rsbmv_test();
    printf("*** Testing Rsbmv successful ***\n");
    return (0);
}
