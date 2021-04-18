/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clarfg.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 2
#define MAX_N 20
#define MIN_INCX  1
#define MAX_INCX  3
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Clarfg_test()
{
    int errorflag = FALSE;
    int incx, n, iter;
    REAL_REF diff;
    COMPLEX_REF alpha_ref, tau_ref;
    COMPLEX_REF *x_ref;
    COMPLEX alpha, tau;
    COMPLEX *x;

    for (incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (n = MIN_N; n <= MAX_N; n++) {
#if defined VERBOSE_TEST
	    printf("# n: %d, incx: %d \n", n, incx);
#endif
	    x_ref = new COMPLEX_REF[veclen(n, incx)];
	    x = new COMPLEX[veclen(n, incx)];
	    for (iter = 0; iter < MAX_ITER; iter++) {
		set_random_vector(x_ref, x, veclen(n, incx));
		set_random_number(alpha_ref, alpha);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		zlarfg_f77(&n, &alpha_ref, x_ref, &incx, &tau_ref);
#else
		Clarfg(n, &alpha_ref, x_ref, incx, &tau_ref);
#endif
		Clarfg(n, &alpha, x, incx, &tau);

		diff = abs(alpha_ref - alpha);
		if (diff > EPSILON) {
		    printf("error alpha: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
	        if (maxdiff < diff)
		    maxdiff = diff;
#if defined VERBOSE_TEST
	        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		diff = abs(tau_ref - tau);
		if (diff > EPSILON) {
		    printf("error tau:   "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
	        if (maxdiff < diff)
		    maxdiff = diff;
#if defined VERBOSE_TEST
	        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
	    }
	    delete[]x_ref;
	    delete[]x;
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Clarfg failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clarfg start ***\n");
    Clarfg_test();
    printf("*** Testing Clarfg successful ***\n");
    return (0);
}
