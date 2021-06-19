/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlaset.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N 10
#define MIN_INC 1
#define MAX_INC 1
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Clacgv_test(void) {
    int errorflag = FALSE;
    int j = 0;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int incx = 1; incx <= MAX_INC; incx++) {
            COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
            COMPLEX *x = new COMPLEX[veclen(n, incx)];
            j = 0;
#if defined VERBOSE_TEST
            printf("#n:%d incx:%d \n", n, incx);
#endif
            while (j < MAX_ITER) {
                set_random_vector(x_ref, x, veclen(n, incx));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                int incx_ref = (int)incx;
                zlacgv_f77(&n, x_ref, &incx_ref);
#else
                Clacgv(n, x_ref, incx);
#endif
                Clacgv(n, x, incx);
                diff = infnorm(x_ref, x, veclen(n, incx), 1);
                if (diff > EPSILON) {
                    errorflag = TRUE;
                    printf("Error\n");
                    printvec(x_ref, n);
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
            delete[] x;
            delete[] x_ref;
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rlaset failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Clacgv start ***\n");
    Clacgv_test();
    printf("*** Testing Clacgv successful ***\n");
    return (0);
}
