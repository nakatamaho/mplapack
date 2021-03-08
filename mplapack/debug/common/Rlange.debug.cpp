/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlange.debug.cpp,v 1.1 2010/08/12 22:52:44 nakatamaho Exp $
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

#define MIN_N     0
#define MAX_N     10
#define MIN_M     0
#define MAX_M     10
#define MAX_LDA   15
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Rlange_test2(const char *norm)
{
    int errorflag = FALSE;
    int j = 0;
    REAL_REF diff;
    REAL_REF Rlange_ref_ret;
    REAL Rlange_ret;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int m = MIN_M; m < MAX_M; m++) {
	    for (int lda = max(1, m); lda < MAX_LDA; lda++) {
		REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
		REAL_REF *work_ref = new REAL_REF[veclen(m, 1)];
		REAL *A = new REAL[matlen(lda, n)];
		REAL *work = new REAL[veclen(m, 1)];

		j = 0;
#if defined VERBOSE_TEST
		printf("#n:%d m:%d lda: %d norm %s\n", n, m, lda, norm);
#endif
		while (j < MAX_ITER) {
		    set_random_vector(A_ref, A, matlen(lda, n));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		    Rlange_ref_ret = dlange_f77(norm, &m, &n, A_ref, &lda, work_ref);
#else
		    Rlange_ref_ret = Rlange(norm, m, n, A_ref, lda, work_ref);
#endif
		    Rlange_ret = Rlange(norm, m, n, A, lda, work);

		    diff = abs(Rlange_ref_ret - Rlange_ret);
		    if (diff > EPSILON) {
			errorflag = TRUE; printf("Error\n"); exit(1);
		    }
	            if (maxdiff < diff) maxdiff = diff;
#if defined VERBOSE_TEST
	            printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		    j++;
		}
		delete[]work_ref;
		delete[]A_ref;
		delete[]work;
		delete[]A;
	    }
	}
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rlange failed ***\n");
	exit(1);
    }
}

void Rlange_test(void)
{
    Rlange_test2("M");
    Rlange_test2("m");
    Rlange_test2("1");
    Rlange_test2("O");
    Rlange_test2("o");
    Rlange_test2("I");
    Rlange_test2("i");
    Rlange_test2("F");
    Rlange_test2("f");
    Rlange_test2("E");
    Rlange_test2("e");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rlange start ***\n");
    Rlange_test();
    printf("*** Testing Rlange successful ***\n");
    return (0);
}
