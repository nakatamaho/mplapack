/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlauum.debug.cpp,v 1.1 2010/08/18 08:50:27 nakatamaho Exp $
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
#include <mplapack.h>
#include <mplapack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N      0
#define MAX_N     25		//should not be so large
#define MAX_LDA   25		//should not be so large
#define MAX_ITER  3

REAL_REF maxdiff = 0.0;

void Rlauum_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
	    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
	    REAL *A = new REAL[matlen(lda, n)];
#if defined VERBOSE_TEST
	    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
#endif
	    j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(A_ref, A, matlen(lda, n));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		dlauum_f77(uplo, &n, A_ref, &lda, &info_ref);
#else
		Rlauum(uplo, n, A_ref, lda, &info_ref);
#endif
		Rlauum(uplo, n, A, lda, &info);

		if (info < 0) {
		    printf("info %d error\n", -(int) info);
		    errorflag = TRUE;
		}
		if (info_ref != info) {
		    printf("info error! %d, %d\n", (int)info_ref, (int)info);
		    errorflag = TRUE;
		}
                diff = infnorm(A_ref, A, matlen(lda, n), 1);
		if (diff > EPSILON10) {
		    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
		    printf("error1: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
	        if (maxdiff < diff) maxdiff = diff;
#if defined VERBOSE_TEST
	        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		j++;
	    }
	    delete[]A;
	    delete[]A_ref;
	}
	if (errorflag == TRUE) {
            printf("*** Testing Rlauum failed ***\n");
	    exit(1);
	}
    }
}

void Rlauum_test(void)
{
    Rlauum_test2("U");
    Rlauum_test2("L");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rlauum start ***\n");
    Rlauum_test();
    printf("*** Testing Rlauum successful ***\n");
    return (0);
}
