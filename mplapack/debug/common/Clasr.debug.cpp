/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clasr.debug.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N     14
#define MIN_M     0
#define MAX_M     14
#define MAX_LDA   16
#define MAX_ITER  2

REAL_REF maxdiff = 0.0;

void Clasr_test2(const char *side, const char *pivot, const char *direct)
{
    int errorflag = FALSE;
    int j = 0;
    int cdim = 0;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int m = MIN_M; m < MAX_M; m++) {
	    for (int lda = max(m, 1); lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
		printf("# n:%d m:%d lda:%d, side %s, pivot %s direct %s\n", n, m, lda, side, pivot, direct);
#endif
		if (Mlsame(side, "L"))
		    cdim = m - 1;
		if (Mlsame(side, "R"))
		    cdim = n - 1;

		COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
		REAL_REF *c_ref = new REAL_REF[max(cdim, 1)];
		REAL_REF *s_ref = new REAL_REF[max(cdim, 1)];
		COMPLEX *A = new COMPLEX[matlen(lda, n)];
		REAL *c = new REAL[max(cdim, 1)];
		REAL *s = new REAL[max(cdim, 1)];

		j = 0;
		while (j < MAX_ITER) {
		    set_random_vector(A_ref, A, matlen(lda, n));
		    set_random_vector(c_ref, c, max(cdim, 1));
		    set_random_vector(s_ref, s, max(cdim, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		    zlasr_f77(side, pivot, direct, &m, &n, c_ref, s_ref, A_ref, &lda);
#else
		    Clasr(side, pivot, direct, m, n, c_ref, s_ref, A_ref, lda);
#endif
		    Clasr(side, pivot, direct, m, n, c, s, A, lda);

                    diff = infnorm(A_ref, A, matlen(lda, n), 1);
		    if (diff > EPSILON) {
			printf("# n:%d m:%d lda:%d, side %s, pivot %s direct %s\n", n, m, lda, side, pivot, direct);
			printf("error: "); printnum(diff); printf("\n");
			errorflag = TRUE;
		    }
	            if (maxdiff < diff)
		        maxdiff = diff;
#if defined VERBOSE_TEST
	            printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		    j++;
		}
		delete[]s_ref;
		delete[]c_ref;
		delete[]A_ref;
		delete[]s;
		delete[]c;
		delete[]A;
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Clasr failed ***\n");
	exit(1);
    }
}

void Clasr_test(void)
{
    Clasr_test2("L", "V", "F");
    Clasr_test2("R", "V", "F");
    Clasr_test2("L", "T", "F");
    Clasr_test2("R", "T", "F");
    Clasr_test2("L", "B", "F");
    Clasr_test2("R", "B", "F");
    Clasr_test2("L", "V", "B");
    Clasr_test2("R", "V", "B");
    Clasr_test2("L", "T", "B");
    Clasr_test2("R", "T", "B");
    Clasr_test2("L", "B", "B");
    Clasr_test2("R", "B", "B");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clasr start ***\n");
    Clasr_test();
    printf("*** Testing Clasr successful ***\n");
    return (0);
}
