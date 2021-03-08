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

#define MIN_INCX  -2
#define MAX_INCX  5
#define MIN_N     0
#define MAX_N     100
#define MAX_ITER 10

void iRamax_test()
{
    int errorflag = FALSE;
    INTEGER_REF idmax;
    mplapackint ifmax;

    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = MIN_N; n < MAX_N; n++) {
#if defined VERBOSE_TEST
	    printf("# n:%d incx:%d \n", n, incx);
#endif
	    REAL_REF *x_ref = new REAL_REF[veclen(n, incx)];
	    REAL *x = new REAL[veclen(n, incx)];
	    int j = 0;

	    while (j < MAX_ITER) {
		set_random_vector(x_ref, x, veclen(n, incx));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
		idmax = idamax_f77(&n, x_ref, &incx);
#else
		idmax = iRamax(n, x_ref, incx);
#endif
	        ifmax = iRamax(n, x, incx);
#if defined VERBOSE_TEST
		cout << "iRamax    :" << ifmax << endl;
		cout << "iRamax_ref:" << idmax << endl;
#endif
		if (ifmax != idmax) {
		    printf("error!!\n");
		    errorflag = TRUE;
		}
		j++;
	    }
	    delete[]x;
	    delete[]x_ref;
	}
    }
    if (errorflag == TRUE) {
	printf("iRamax test failed...\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing iRamax start ***\n");
    iRamax_test();
    printf("*** Testing iRamax successful ***\n");
    return (0);
}
