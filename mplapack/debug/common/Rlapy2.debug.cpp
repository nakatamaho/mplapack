/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlapy2.debug.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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

REAL_REF maxdiff = 0.0;

void Rlapy2_test()
{
    int errorflag = FALSE;
    REAL_REF dlapy2_ret, x_ref, y_ref, diff;
    REAL Rlapy2_ret, x, y;

    int count = 100;

    while (count--) {
	set_random_number(x_ref, x);
	set_random_number(y_ref, y);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
	dlapy2_ret = dlapy2_f77(&x_ref, &y_ref);
#else
	dlapy2_ret = Rlapy2(x_ref, y_ref);
#endif
	Rlapy2_ret = Rlapy2(x, y);

#if defined VERBOSE_TEST
	cout << "x      "; printnum(x); cout << endl;
	cout << "y      "; printnum(y); cout << endl;
	cout << "Rlapy2 "; printnum(Rlapy2_ret); cout << endl;
	cout << "x_ref  "; printnum(x_ref); cout << endl;
	cout << "y_ref  "; printnum(y_ref); cout << endl;
	cout << "dlapy2 "; printnum(dlapy2_ret); cout << endl;
#endif
        diff = abs (dlapy2_ret - Rlapy2_ret);
	if (diff > EPSILON) {
	    errorflag = TRUE;
	}
	if (maxdiff < diff)
	    maxdiff = diff;
#if defined VERBOSE_TEST
	printf("max error: "); printnum(maxdiff); printf("\n");
#endif
    }
    if (errorflag == TRUE) {
	printf("*** Testing Rlapy2 failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rlapy2 start ***\n");
    Rlapy2_test();
    printf("*** Testing Rlapy2 successful ***\n");
    return (0);
}
