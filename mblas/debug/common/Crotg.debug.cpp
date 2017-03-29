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
#include <mblas.h>
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MAX_ITER 1000

REAL_REF maxdiff = 0.0;

void Crotg_test()
{
    int errorflag = FALSE;

    COMPLEX_REF ca_ref;
    COMPLEX_REF cb_ref;
    REAL_REF cc_ref;
    COMPLEX_REF cs_ref;
    REAL_REF diff1;
    REAL_REF diff2;
    REAL_REF diff3;
    REAL_REF diff4;

    COMPLEX ca;
    COMPLEX cb;
    REAL cc;
    COMPLEX cs;

    int j = 0;
    while (j < MAX_ITER) {
	set_random_number(ca_ref, ca);
	set_random_number(cb_ref, cb);
	set_random_number(cc_ref, cc);

#if defined ___MPACK_BUILD_WITH_MPFR___
	zrotg_f77(&ca_ref, &cb_ref, &cc_ref, &cs_ref);
#else
	Crotg(&ca_ref, cb_ref, &cc_ref, &cs_ref);
#endif
	Crotg(&ca, cb, &cc, &cs);

	diff1 = abs(ca_ref - ca);
	diff2 = abs(cb_ref - cb);
	diff3 = abs(cc_ref - cc);
	diff4 = abs(cs_ref - cs);

#if defined VERBOSE_TEST
	printf("diff1="); printnum(diff1); printf("\n");
	printf("diff2="); printnum(diff2); printf("\n");
	printf("diff3="); printnum(diff3); printf("\n");
	printf("diff4="); printnum(diff4);printf("\n");
#endif
	if (maxdiff < diff1) maxdiff = diff1;
	if (maxdiff < diff2) maxdiff = diff2;
	if (maxdiff < diff3) maxdiff = diff3;
	if (maxdiff < diff4) maxdiff = diff4;

	if (diff1 > EPSILON) {
#if defined VERBOSE_TEST
	    printf("error: "); printnum(diff1); printf("\n");
#endif
	    errorflag = TRUE;
	}
	if (diff2 > EPSILON) {
#if defined VERBOSE_TEST
	    printf("error: "); printnum(diff2); printf("\n");
#endif
	    errorflag = TRUE;
	}
	if (diff3 > EPSILON) {
#if defined VERBOSE_TEST
	    printf("error: "); printnum(diff3); printf("\n");
#endif
	    errorflag = TRUE;
	}
	if (diff4 > EPSILON) {
#if defined VERBOSE_TEST
	    printf("error: "); printnum(diff4); printf("\n");
#endif
	    errorflag = TRUE;
	}
	j++;
    }

    if (errorflag == TRUE) {
	printf("error: "); printnum(diff1); printf("\n");
	printf("error: "); printnum(diff2); printf("\n");
	printf("error: "); printnum(diff3); printf("\n");
	printf("error: "); printnum(diff4); printf("\n");
	printf("*** Testing Crotg failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Crotg start ***\n");
    Crotg_test();
    printf("*** Testing Crotg successful ***\n");
    return (0);
}
