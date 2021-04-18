/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clatrd.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_NB 1
#define MAX_NB 10
#define MIN_N 4
#define MAX_N 10
#define MIN_LDA 4
#define MAX_LDA 10
#define MIN_LDW 4
#define MAX_LDW 10
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Clatrd_test2(const char *uplo)
{
    int errorflag = FALSE;
    int iter;
    int n, nb;
    INTEGER_REF lda, ldw;
    REAL_REF diff;

    for (n = MIN_N; n <= MAX_N; n++) {
      for (lda = max(1, n); lda <= MAX_LDA; lda++) {
	  for (nb = MIN_NB; nb <= lda && nb <=n ; nb++) {
		for (ldw = max(1, n); ldw <= MAX_LDW; ldw++) {
#if defined VERBOSE_TEST
		    printf("# uplo %s, n %d, lda %d, nb %d\n", uplo, n, (int)lda, nb);
#endif
		    COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
		    COMPLEX_REF *W_ref = new COMPLEX_REF[matlen(ldw, nb)];
		    REAL_REF *e_ref = new REAL_REF[veclen(n - 1, 1)];
		    COMPLEX_REF *tau_ref = new COMPLEX_REF[veclen(n - 1, 1)];

		    COMPLEX *A = new COMPLEX[matlen(lda, n)];
		    COMPLEX *W = new COMPLEX[matlen(ldw, nb)];
		    REAL *e = new REAL[veclen(n - 1, 1)];
		    COMPLEX *tau = new COMPLEX[veclen(n - 1, 1)];

		    for (iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(A_ref, A, matlen(lda, n));
			set_random_vector(W_ref, W, matlen(ldw, nb));
			set_random_vector(e_ref, e, veclen(n - 1, 1));
			set_random_vector(tau_ref, tau, veclen(n - 1, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
			zlatrd_f77(uplo, &n, &nb, A_ref, &lda, e_ref, tau_ref, W_ref, &ldw);
#else
			Clatrd(uplo, n, nb, A_ref, lda, e_ref, tau_ref, W_ref, ldw);
#endif
			Clatrd(uplo, n, nb, A, lda, e, tau, W, ldw);

			diff = infnorm(A_ref, A, matlen(lda, n), 1);
		        if (diff > EPSILON12) {
		            printf("error in A: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;

			diff = infnorm(W_ref, W, matlen(ldw, nb), 1);
		        if (diff > EPSILON12) {
		            printf("error in W: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;

			diff = infnorm(e_ref, e, veclen(n - 1, 1), 1);
		        if (diff > EPSILON12) {
		            printf("error in e: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;

			diff = infnorm(tau_ref, tau, veclen(n - 1, 1), 1);
		        if (diff > EPSILON12) {
		            printf("error in tau:"); printnum(diff); printf("\n");
		            errorflag = TRUE;
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;
#if defined VERBOSE_TEST
                        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		    }
		    delete[]tau_ref;
		    delete[]e_ref;
		    delete[]W_ref;
		    delete[]A_ref;
		    delete[]tau;
		    delete[]e;
		    delete[]W;
		    delete[]A;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Clatrd failed ***\n");
	exit(1);
    }
}

void Clatrd_test()
{
    Clatrd_test2("L");
    Clatrd_test2("U");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clatrd start ***\n");
    Clatrd_test();
    printf("*** Testing Clatrd successful ***\n");
    return (0);
}
