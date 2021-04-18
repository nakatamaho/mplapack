/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 * $Id: Rlasrt.debug.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N     50
#define MAX_LDA   60
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Rlasrt_test2(const char *id)
{
    int errorflag = FALSE;
    int j;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
#if defined VERBOSE_TEST
	printf("#n:%d id %s\n", n, id);
#endif
	REAL_REF *d_ref = new REAL_REF[veclen(n, 1)];
	REAL *d = new REAL[veclen(n, 1)];

	j = 0;
	while (j < MAX_ITER) {
	    set_random_vector(d_ref, d, veclen(n, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
	    dlasrt_f77(id, &n, d_ref, &info_ref);
#else
	    Rlasrt(id, n, d_ref, &info_ref);
#endif
	    Rlasrt(id, n, d, &info);

	    if (info_ref != info) {
		printf("#error info info:%d info:%d\n", (int) info_ref, (int)info);
		errorflag = TRUE;
	    }
	    diff = infnorm(d_ref, d, veclen(n, 1), 1);
	    if (diff > EPSILON) {
	        printf("error2: "); printnum(diff); printf("\n");
	        errorflag = TRUE;
                  exit(1);
	    }
	    if (maxdiff < diff)
	        maxdiff = diff;
#if defined VERBOSE_TEST
	    printf("max error: "); printnum(maxdiff); printf("\n");
#endif
	    j++;
        }
	delete[]d;
	delete[]d_ref;
    }
    if (errorflag == TRUE) {
	printf("*** Testing Rlasrt start ***\n");
	exit(1);
    }
}

void Rlasrt_test(void)
{
    Rlasrt_test2("I");
    Rlasrt_test2("D");
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rlasrt start ***\n");
    Rlasrt_test();
    printf("*** Testing Rlasrt successful ***\n");
    return (0);
}
