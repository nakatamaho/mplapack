/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clarft.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_K 1
#define MAX_K 12
#define MIN_N 0
#define MAX_N 12
#define MAX_LDT 1
#define MAX_LDC 12
#define MAX_LDV 12
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

void Clarft_test2(const char *direct, const char *storev)
{
    int errorflag = FALSE;
    int n, k, ldvmin, v_column, ldt, iter, ldv;
    REAL_REF diff;

    for (k = MIN_K; k <= MAX_K; k++) {
	for (n = k; n <= MAX_N; n++) {
	    if (Mlsame(storev, "C"))
		ldvmin = max(1, n);
	    else //storev = "R"
		ldvmin = k;
	    for (ldv = ldvmin; ldv <= MAX_LDV; ldv++) {
		for (ldt = k; ldt <= MAX_LDT; ldt++) {
		    if (Mlsame(storev, "C"))
			v_column = k;
		    else //storev = "R"
			v_column = n;
#if defined VERBOSE_TEST
		    printf("#direct %s: storev %s, n %d, k %d, ldv %d, ldt %d\n", direct, storev, n, k, ldv, ldt);

#endif
		    COMPLEX_REF *T_ref = new COMPLEX_REF[matlen(ldt, k)];
		    COMPLEX_REF *V_ref = new COMPLEX_REF[matlen(ldv, v_column)];
		    COMPLEX_REF *tau_ref = new COMPLEX_REF[veclen(k, 1)];

		    COMPLEX *T = new COMPLEX[matlen(ldt, k)];
		    COMPLEX *V = new COMPLEX[matlen(ldv, v_column)];
		    COMPLEX *tau = new COMPLEX[veclen(k, 1)];

		    for (iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(T_ref, T, matlen(ldt, k));
			set_random_vector(V_ref, V, matlen(ldv, v_column));
			set_random_vector(tau_ref, tau, veclen(k, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			zlarft_f77(direct, storev, &n, &k, V_ref, &ldv, tau_ref, T_ref, &ldt);
#else
			Clarft(direct, storev, n, k, V_ref, ldv, tau_ref, T_ref, ldt);
#endif
			Clarft(direct, storev, n, k, V, ldv, tau, T, ldt);

			diff = infnorm(T_ref, T, matlen(ldt, k), 1);
		        if (diff > EPSILON) {
		            printf("error in T: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;

			diff = infnorm(V_ref, V, matlen(ldv, v_column), 1);
		        if (diff > EPSILON) {
		            printf("error in V: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;
#if defined VERBOSE_TEST
                        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		    }
		    delete[]tau_ref;
		    delete[]V_ref;
		    delete[]T_ref;
		    delete[]tau;
		    delete[]V;
		    delete[]T;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Clarft failed ***\n");
	exit(1);
    }
}

void Clarft_test()
{
    Clarft_test2("B", "C");
    Clarft_test2("B", "R");
    Clarft_test2("F", "C");
    Clarft_test2("F", "R");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clarft start ***\n");
    Clarft_test();
    printf("*** Testing Clarft successful ***\n");
    return (0);
}
