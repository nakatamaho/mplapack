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

#define MIN_K -2
#define MAX_K 6
#define MIN_N -2
#define MAX_N 6
#define MIN_LDA -2
#define MAX_LDA 6
#define MIN_LDB -2
#define MAX_LDB 6
#define MIN_LDC -2
#define MAX_LDC 6
#define MAX_ITER 1

REAL_REF maxdiff = 0.0;

void Cher2k_test3(const char *uplo, const char *trans, COMPLEX_REF alpha_ref, REAL_REF beta_ref, COMPLEX alpha, REAL beta)
{
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int n = MIN_N; n < MAX_N; n++) {
	for (int k = MIN_K; k < MAX_K; k++) {
	    int minlda;
	    if (Mlsame(trans, "N"))
		minlda = max(1, n);
	    else
		minlda = max(1, k);
	    int minldb;
	    if (Mlsame(trans, "N"))
		minldb = max(1, n);
	    else
		minldb = max(1, k);
	    for (int lda = minlda; lda < MAX_LDA; lda++) {
		for (int ldb = minldb; ldb < MAX_LDB; ldb++) {
		    for (int ldc = max(1, n); ldc < MAX_LDC; ldc++) {
#if defined VERBOSE_TEST
			printf("#n is %d, k is %d, lda is %d, ldb is %d, ldc is %d ", n, k, lda, ldb, ldc);
			printf("uplo is %s, trans is %s \n", uplo, trans);
#endif
			COMPLEX_REF *A_ref;
			COMPLEX_REF *B_ref;
			COMPLEX_REF *C_ref;
			COMPLEX *A;
			COMPLEX *B;
			COMPLEX *C;
			int ka, kb;
			if (Mlsame(trans, "N"))
			    ka = k;
			else
			    ka = n;
			if (Mlsame(trans, "N"))
			    kb = k;
			else
			    kb = n;

			A_ref = new COMPLEX_REF[matlen(lda, ka)];
			B_ref = new COMPLEX_REF[matlen(ldb, kb)];
			C_ref = new COMPLEX_REF[matlen(ldc, n)];
			A = new COMPLEX[matlen(lda, ka)];
			B = new COMPLEX[matlen(ldb, kb)];
			C = new COMPLEX[matlen(ldc, n)];

			for (int iter = 0; iter < MAX_ITER; iter++) {
			    set_random_vector(A_ref, A, matlen(lda, ka));
			    set_random_vector(B_ref, B, matlen(ldb, kb));
			    set_random_vector(C_ref, C, matlen(ldc, n));

			    mplapack_errno = 0; blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			    zher2k_f77(uplo, trans, &n, &k, &alpha_ref, A_ref, &lda, B_ref, &ldb, &beta_ref, C_ref, &ldc);
			    mplapack_errno1 = blas_errno;
#else
			    Cher2k(uplo, trans, n, k, alpha_ref, A_ref, lda, B_ref, ldb, beta_ref, C_ref, ldc);
			    mplapack_errno1 = mplapack_errno;
#endif
			    Cher2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
			    mplapack_errno2 = mplapack_errno;
#if defined VERBOSE_TEST
			    printf("errno: mplapack %d, ref %d\n", mplapack_errno1,  mplapack_errno2);
#endif
			    if (mplapack_errno1 != mplapack_errno2) {
#if defined VERBOSE_TEST
				printf("error in Mxerbla!!\n");
#endif
				errorflag = TRUE;
			    }
			    REAL_REF diff = infnorm(C_ref, C, matlen(ldc, n), 1);
			    if (diff > EPSILON) {
#if defined VERBOSE_TEST
				printf("error: "); printnum(diff); printf("\n");
#endif
				errorflag = TRUE;
			    }
			    if (maxdiff < diff)
				maxdiff = diff;
			    delete[]C_ref;
			    delete[]B_ref;
			    delete[]A_ref;
			    delete[]C;
			    delete[]B;
			    delete[]A;
			}
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Cher2k failed ***\n");
	exit(1);
    }
}

void Cher2k_test2(const char *uplo, const char *trans)
{
    COMPLEX_REF alpha_ref;
    REAL_REF beta_ref;
    COMPLEX alpha;
    REAL beta;

//a=0, b=*;
    alpha_ref = 0.0;
    alpha = 0.0;
    set_random_number(beta_ref, beta);
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);
//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0;
    beta = 0.0;
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0;
    beta = 1.0;
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0; beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0;
    alpha = 1.0;
    set_random_number(beta_ref, beta);
    Cher2k_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);
}

void Cher2k_test()
{
    Cher2k_test2("U", "N");
    Cher2k_test2("U", "C");
    Cher2k_test2("L", "N");
    Cher2k_test2("L", "C");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Cher2k start ***\n");
    Cher2k_test();
    printf("*** Testing Cher2k successful ***\n");
    return (0);
}
