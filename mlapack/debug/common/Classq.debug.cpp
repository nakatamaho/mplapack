/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Classq.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_INCX  1
#define MAX_INCX  10
#define MIN_N     1
#define MAX_N     100
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Classq_test()
{
    int errorflag = FALSE;
    REAL_REF scale_ref, sumsq_ref, diff;
    REAL scale, sumsq;

    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = MIN_N; n <= MAX_N; n++) {
	    COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
	    COMPLEX *x = new COMPLEX[veclen(n, incx)];
#if defined VERBOSE_TEST
	    printf("# incx: %d, n: %d\n", incx, n);
#endif
	    int j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(x_ref, x, veclen(n, incx));
		set_random_number(scale_ref, scale);
		set_random_number(sumsq_ref, sumsq);
#if defined ___MPACK_BUILD_WITH_MPFR___
		zlassq_f77(&n, x_ref, &incx, &scale_ref, &sumsq_ref);
#else
		Classq(n, x_ref, incx, &scale_ref, &sumsq_ref);
#endif
		Classq(n, x, incx, &scale, &sumsq);

		diff = abs(scale_ref - scale);
		if (diff > EPSILON) {
 		    printf("error: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
                if (maxdiff < diff)
                    maxdiff = diff;
#if defined VERBOSE_TEST
 	        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		diff = abs(sumsq_ref - sumsq);
		if (diff > EPSILON) {
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
	    delete[]x;
	    delete[]x_ref;
	}
    }
    if (errorflag == TRUE) {
        printf("*** Testing Classq failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Classq start ***\n");
    Classq_test();
    printf("*** Testing Classq successful ***\n");
    return (0);
}
