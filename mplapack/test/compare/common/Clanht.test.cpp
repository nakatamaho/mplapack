/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clanht.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
#include <mplapack_compare_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N 0
#define MAX_N 100
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Clanht_test2(const char *norm) {
    int errorflag = FALSE;
    int j = 0, ni;
    REAL_REF zlanht_ret, diff;
    REAL Clanht_ret;

    for (mplapackint n = MIN_N; n < MAX_N; n++) {
        ni = n;
#if defined VERBOSE_TEST
        printf("n:%d norm %s\n", (int)n, norm);
#endif
        REAL_REF *d_ref = new REAL_REF[veclen(n, 1)];
        COMPLEX_REF *e_ref = new COMPLEX_REF[veclen(n - 1, 1)];
        REAL *d = new REAL[veclen(n, 1)];
        COMPLEX *e = new COMPLEX[veclen(n - 1, 1)];
        j = 0;
        while (j < MAX_ITER) {
            set_random_vector(d_ref, d, veclen(n, 1));
            set_random_vector(e_ref, e, veclen(n - 1, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
            zlanht_ret = zlanht_f77(norm, &ni, d_ref, e_ref);
#else
            zlanht_ret = Clanht(norm, n, d_ref, e_ref);
#endif
            Clanht_ret = Clanht(norm, n, d, e);

            diff = abs(Clanht_ret - zlanht_ret);
            if (diff > EPSILON) {
                printf("error: ");
                printnum(diff);
                printf("\n");
                errorflag = TRUE;
            }
            if (maxdiff < diff)
                maxdiff = diff;
#if defined VERBOSE_TEST
            printf("max error: ");
            printnum(maxdiff);
            printf("\n");
#endif
            j++;
        }
        delete[] d;
        delete[] d_ref;
        delete[] e;
        delete[] e_ref;
    }
    if (errorflag == TRUE) {
        printf("*** Testing Clanht failed ***\n");
        exit(1);
    }
#if defined VERBOSE_TEST
    printf("max error: ");
    printnum(maxdiff);
    printf("\n");
#endif
}

void Clanht_test(void) {
    Clanht_test2("M");
    Clanht_test2("m");
    Clanht_test2("1");
    Clanht_test2("O");
    Clanht_test2("o");
    Clanht_test2("I");
    Clanht_test2("i");
    Clanht_test2("F");
    Clanht_test2("f");
    Clanht_test2("E");
    Clanht_test2("e");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Clanht start ***\n");
    Clanht_test();
    printf("*** Testing Clanht successful ***\n");
    return (0);
}
