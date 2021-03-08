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

#define MAX_ITER  100

REAL_REF maxdiff = 0.0;

void Rrotg_test()
{
    int errorflag = FALSE;
    REAL_REF a_ref;
    REAL_REF b_ref;
    REAL_REF c_ref;
    REAL_REF s_ref;
    REAL_REF diff[4];
    REAL a;
    REAL b;
    REAL c;
    REAL s;
    int j = 0;
    while (j < MAX_ITER) {
	set_random_number(a_ref, a);
	set_random_number(b_ref, b);
	set_random_number(c_ref, c);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
	drotg_f77(&a_ref, &b_ref, &c_ref, &s_ref);
#else
	Rrotg(&a_ref, &b_ref, &c_ref, &s_ref);
#endif
	Rrotg(&a, &b, &c, &s);
	diff[0] = abs(a - a_ref);
	diff[1] = abs(b - b_ref);
	diff[2] = abs(c - c_ref);
	diff[3] = abs(s - s_ref);
#if defined VERBOSE_TEST
	for (int p = 0; p < 4; p++) {
	    printf("diff[%d]=", p); printnum(diff[p]); printf("\n");
	}
#endif
	for (int p = 0; p < 4; p++) {
	    if (maxdiff < diff[p])
	   	 maxdiff = diff[p];
	    if (diff[p] > EPSILON2) {
#if defined VERBOSE_TEST
		printf("error %d,\n", p); printnum(diff[p]); printf("\n");
#endif
		errorflag = TRUE;
	    }
	}
	j++;
    }
    if (errorflag == TRUE) {
        printf("error: "); printnum(maxdiff); printf("\n");
        printf("*** Testing Rrotg failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rrotg start ***\n");
    Rrotg_test();
    printf("*** Testing Rrotg successful ***\n");
    return (0);
}
