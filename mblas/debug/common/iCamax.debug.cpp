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

#define MIN_INCX  1
#define MAX_INCX  2
#define MIN_N     1
#define MAX_N     10
#define MAX_ITER 10

void iCamax_test()
{
    int errorflag = FALSE;
    INTEGER_REF izmax;
    INTEGER icmax;

    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = MIN_N; n < MAX_N; n++) {
#if defined VERBOSE_TEST
	    printf("#n:%d incx:%d \n", n, incx);
#endif
	    COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
	    COMPLEX *x = new COMPLEX[veclen(n, incx)];
	    int j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(x_ref, x, veclen(n, incx));
#if defined ___MPACK_BUILD_WITH_MPFR___
		izmax = izamax_f77(&n, x_ref, &incx);
#else
		izmax = iCamax(n, x_ref, incx);
#endif
		icmax = iCamax(n, x, incx);
#if defined VERBOSE_TEST
		cout << "iCamax:" << icmax << endl;
		cout << "izamax:" << izmax << endl;
#endif
		if (icmax != izmax) {
#if defined VERBOSE_TEST
		    printf("icmax!=icmax_ref\n");
#endif
		    errorflag = TRUE;
		}
		j++;
	    }
	    delete[]x;
	    delete[]x_ref;
	}
    }
    if (errorflag == TRUE) {
	printf("*** Testing iCamax failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing iCamax start ***\n");
    iCamax_test();
    printf("*** Testing iCamax successful ***\n");
    return (0);
}
