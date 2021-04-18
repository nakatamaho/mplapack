/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clacrt.debug.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MAX_INCX  10
#define MIN_INCX -10
#define MAX_INCY  10
#define MIN_INCY -10
#define MAX_N     10
#define MAX_ITER  1

REAL_REF maxdiff = 0.0;

void Clacrt_test()
{
    int errorflag = FALSE;
    REAL_REF diff;
    COMPLEX ctemp;
    for (int incx = MAX_INCX; incx >= MIN_INCX; incx--) {
	for (int incy = MAX_INCY; incy >= MIN_INCY; incy--) {
	    for (int n = 3; n < MAX_N; n++) {
#if defined VERBOSE
		printf("# n:%d incx:%d, incy:%d\n", n, incx, incy);
#endif
		COMPLEX_REF *cx_ref = new COMPLEX_REF[veclen(n, incx)];
		COMPLEX_REF *cy_ref = new COMPLEX_REF[veclen(n, incy)];
		COMPLEX_REF c_ref, s_ref;

		COMPLEX *cx = new COMPLEX[veclen(n, incx)];
		COMPLEX *cy = new COMPLEX[veclen(n, incy)];
		COMPLEX c, s;

		set_random_number(c_ref, c);
		set_random_number(s_ref, s);
		int j = 0;
		while (j < MAX_ITER) {
		    set_random_vector(cx_ref, cx, veclen(n, incx));
		    set_random_vector(cy_ref, cy, veclen(n, incy));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		    zlacrt_f77(&n, cx_ref, &incx, cy_ref, &incy, &c_ref, &s_ref);
#else
		    Clacrt(n, cx_ref, incx, cy_ref, incy, c_ref, s_ref);
#endif
		    Clacrt(n, cx, incx, cy, incy, c, s);

		    diff = infnorm(cx_ref, cx, veclen(n, incx), 1);
		    if (diff > EPSILON) {
			printf("error: "); printnum(diff); printf("\n");
			errorflag = TRUE;
		    }
		    if (maxdiff < diff)
			maxdiff = diff;
#if defined VERBOSE
		    printf("max error: "); printnum(maxdiff); printf("\n");
#endif
		    j++;
		}
		delete[]cx;
		delete[]cy;
		delete[]cx_ref;
		delete[]cy_ref;
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing Clacrt failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Clacrt start ***\n");
    Clacrt_test();
    printf("*** Testing Clacrt successful ***\n");
    return (0);
}
