/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rladiv.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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

REAL_REF maxdiff = 0.0;

void Rladiv_test()
{
    int errorflag = FALSE;
    REAL_REF a_ref, b_ref, c_ref, d_ref, p_ref, q_ref, diff;
    REAL a, b, c, d, p, q;

    int count = 100;
    while (count--) {
	set_random_number(a_ref, a);
	set_random_number(b_ref, b);
	set_random_number(c_ref, c);
	set_random_number(d_ref, d);

#if defined ___MPACK_BUILD_WITH_MPFR___
	dladiv_f77(&a_ref, &b_ref, &c_ref, &d_ref, &p_ref, &q_ref);
#else
	Rladiv(a_ref, b_ref, c_ref, d_ref, &p_ref, &q_ref);
#endif
        Rladiv(a, b, c, d, &p, &q);

#if defined VERBOSE_TEST
	cout << "a      "; printnum(a); cout << endl;
	cout << "b      "; printnum(b); cout << endl;
	cout << "c      "; printnum(c); cout << endl;
	cout << "d      "; printnum(d); cout << endl;
	cout << "p      "; printnum(p);	cout << endl;
	cout << "q      "; printnum(q); cout << endl;

	cout << "a_ref  "; printnum(a_ref); cout << endl;
	cout << "b_ref  "; printnum(b_ref); cout << endl;
	cout << "c_ref  "; printnum(c_ref); cout << endl;
	cout << "d_ref  "; printnum(d_ref); cout << endl;
	cout << "p_ref  "; printnum(p_ref); cout << endl;
	cout << "q_ref  "; printnum(q_ref); cout << endl;
#endif
        diff = abs (p_ref - p);
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("error1: "); printnum(diff); printf("\n");
	}
        if (maxdiff < diff)
	    maxdiff = diff;
        diff = abs (q_ref - q);
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("error2: "); printnum(diff); printf("\n");
	}
        if (maxdiff < diff)
	    maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: "); printnum(maxdiff); printf("\n");
#endif
    }
    if (errorflag == TRUE) {
	printf("*** Testing Rladiv failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Rladiv start ***\n");
    Rladiv_test();
    printf("*** Testing Rladiv successful ***\n");
    return (0);
}
