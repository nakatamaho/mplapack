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

#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void RCabs1_test()
{
    int errorflag = FALSE;
    COMPLEX_REF z_ref;
    REAL_REF dtemp;

    COMPLEX z;
    REAL Ftemp;

    int j = 0;

    while (j < MAX_ITER) {
	set_random_number(z_ref, z);

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
	dtemp = dcabs1_f77(&z_ref);
#else
	dtemp = RCabs1(z_ref);
#endif
	Ftemp = RCabs1(z);

#if defined VERBOSE_TEST
	printf("Z     ="); printnum(z); printf("\n");
	printf("z_ref ="); printnum(z_ref); printf("\n");

	printf("RCabs1 = "); printnum(Ftemp); printf("\n");
	printf("dcabs1 = "); printnum(dtemp); printf("\n");
#endif
	REAL_REF diff = abs(Ftemp - dtemp);
	if (diff > EPSILON) {
#if defined VERBOSE_TEST
	    printf("error: "); printnum(diff); printf("\n");
#endif
	    errorflag = TRUE;
	}
	if (maxdiff < diff)
	    maxdiff = diff;
	j++;
    }
    if (errorflag == TRUE) {
	printf("error: "); printnum(maxdiff); printf("\n");
	printf("*** Testing RCabs1 failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing RCabs1 start ***\n");
    RCabs1_test();
    printf("*** Testing RCabs1 successful ***\n");
    return (0);
}
