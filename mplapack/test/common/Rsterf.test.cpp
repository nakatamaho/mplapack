/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rsterf.debug.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 0
#define MAX_N 20
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rsterf_test(void) {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    REAL_REF diff;
    INTEGER info;

    for (int n = MIN_N; n < MAX_N; n++) {
        REAL_REF *D_ref = new REAL_REF[veclen(n, 1)];
        REAL_REF *E_ref = new REAL_REF[veclen(n - 1, 1)];
        REAL *D = new REAL[veclen(n, 1)];
        REAL *E = new REAL[veclen(n - 1, 1)];
#if defined VERBOSE_TEST
        printf("# n:%d\n", n);
#endif
        j = 0;
        while (j < MAX_ITER) {
            set_random_vector(D_ref, D, veclen(n, 1));
            set_random_vector(E_ref, E, veclen(n - 1, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
            dsterf_f77(&n, D_ref, E_ref, &info_ref);
#else
            Rsterf(n, D_ref, E_ref, &info_ref);
#endif
            Rsterf(n, D, E, &info);

            if (info < 0) {
                printf("info %d error\n", -(int)info);
                errorflag = TRUE;
            }
            if (info_ref != info) {
                printf("info differ! %d, %d\n", (int)info_ref, (int)info);
                errorflag = TRUE;
            }
            diff = infnorm(D_ref, D, veclen(n, 1), 1);
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
        delete[] D;
        delete[] E;
        delete[] D_ref;
        delete[] E_ref;
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rsterf failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rsterf start ***\n");
    Rsterf_test();
    printf("*** Testing Rsterf successful ***\n");
    return (0);
}
