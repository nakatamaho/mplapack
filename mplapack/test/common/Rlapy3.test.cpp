/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlapy3.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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

void Rlapy3_test() {
    int errorflag = FALSE;
    REAL_REF dlapy3_ret, x_ref, y_ref, z_ref, diff;
    REAL Rlapy3_ret, x, y, z;

    int count = 100;

    while (count--) {
        set_random_number(x_ref, x);
        set_random_number(y_ref, y);
        set_random_number(z_ref, z);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
        dlapy3_ret = dlapy3_f77(&x_ref, &y_ref, &z_ref);
#else
        dlapy3_ret = Rlapy3(x_ref, y_ref, z_ref);
#endif
        Rlapy3_ret = Rlapy3(x, y, z);

#if defined VERBOSE_TEST
        cout << "x      ";
        printnum(x);
        cout << endl;
        cout << "y      ";
        printnum(y);
        cout << endl;
        cout << "z      ";
        printnum(z);
        cout << endl;
        cout << "Rlapy3 ";
        printnum(Rlapy3_ret);
        cout << endl;
        cout << "x_ref  ";
        printnum(x_ref);
        cout << endl;
        cout << "y_ref  ";
        printnum(y_ref);
        cout << endl;
        cout << "z_ref  ";
        printnum(z_ref);
        cout << endl;
        cout << "dlapy3 ";
        printnum(dlapy3_ret);
        cout << endl;
#endif
        diff = abs(dlapy3_ret - Rlapy3_ret);
        if (diff > EPSILON) {
            errorflag = TRUE;
            printf("Error1\n");
        }
        if (maxdiff < diff)
            maxdiff = diff;
#if defined VERBOSE_TEST
        printf("max error: ");
        printnum(maxdiff);
        printf("\n");
#endif
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rlapy3 failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rlapy3 start ***\n");
    Rlapy3_test();
    printf("*** Testing Rlapy3 successful ***\n");
    return (0);
}
