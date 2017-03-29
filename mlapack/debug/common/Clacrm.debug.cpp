/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clacrm.debug.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 1
#define MAX_N 10
#define MIN_M 1
#define MAX_M 10
#define MAX_LDA 10
#define MAX_LDB 10
#define MAX_LDC 10
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Clacrm_test()
{
    int errorflag = FALSE;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int m = MIN_M; m < MAX_M; m++) {
	    for (int lda = max(1, m); lda < MAX_LDA; lda++) {
		for (int ldb = max(1, n); ldb < MAX_LDB; ldb++) {
		    for (int ldc = max(1, m); ldc < MAX_LDC; ldc++) {
#if defined VERBOSE_TEST
			printf("# n %d m %d lda %d ldb %d ldc %d\n", n, m, lda, ldb, ldc);
#endif
			COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
			REAL_REF *B_ref = new REAL_REF[matlen(ldb, n)];
			COMPLEX_REF *C_ref = new COMPLEX_REF[matlen(ldc, n)];
			REAL_REF *rwork_ref = new REAL_REF[veclen(2 * m * n, 1)];

			COMPLEX *A = new COMPLEX[matlen(lda, n)];
			REAL *B = new REAL[matlen(ldb, n)];
			COMPLEX *C = new COMPLEX[matlen(ldc, n)];
			REAL *rwork = new REAL[veclen(2 * m * n, 1)];

			for (int iter = 0; iter < MAX_ITER; iter++) {
			    set_random_vector(A_ref, A, matlen(lda, n));
			    set_random_vector(B_ref, B, matlen(ldb, n));
			    set_random_vector(C_ref, C, matlen(ldc, n));
#if defined ___MPACK_BUILD_WITH_MPFR___
			    zlacrm_f77(&m, &n, A_ref, &lda, B_ref, &ldb, C_ref, &ldc, rwork_ref);
#else
			    Clacrm(m, n, A_ref, lda, B_ref, ldb, C_ref, ldc, rwork_ref);
#endif
			    Clacrm(m, n, A, lda, B, ldb, C, ldc, rwork);

			    diff = infnorm(C_ref, C, matlen(ldc, n), 1);
			    if (diff > EPSILON) {
				printf("error: "); printnum(diff); printf("\n");
				errorflag = TRUE;
			    }
			    if (maxdiff < diff)
				maxdiff = diff;
#if defined VERBOSE_TEST
			    printf("max error: "); printnum(maxdiff); printf("\n");
#endif
			}
			delete[]rwork_ref;
			delete[]C_ref;
			delete[]B_ref;
			delete[]A_ref;
			delete[]rwork;
			delete[]C;
			delete[]B;
			delete[]A;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Clacrm failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clacrm start ***\n");
    Clacrm_test();
    printf("*** Testing Clacrm successful ***\n");
    return (0);
}
