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
#include <mplapack_compare_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MAX_INCX 10
#define MIN_INCX -10
#define MAX_INCY 10
#define MIN_INCY -10
#define MAX_N 10
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Cdotu_test() {
    int errorflag = FALSE;
    COMPLEX_REF dtemp;
    COMPLEX ctemp;
    for (int incx = MAX_INCX; incx >= MIN_INCX; incx--) {
        for (int incy = MAX_INCY; incy >= MIN_INCY; incy--) {
            for (int n = 0; n < MAX_N; n++) {
#if defined VERBOSE_TEST
                printf("# n:%d incx:%d, incy:%d\n", n, incx, incy);
#endif
                COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
                COMPLEX_REF *y_ref = new COMPLEX_REF[veclen(n, incy)];
                COMPLEX *x = new COMPLEX[veclen(n, incx)];
                COMPLEX *y = new COMPLEX[veclen(n, incy)];
                int j = 0;
                while (j < MAX_ITER) {
                    set_random_vector(x_ref, x, veclen(n, incx));
                    set_random_vector(y_ref, y, veclen(n, incy));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    dtemp = zdotu_f77(&n, x_ref, &incx, y_ref, &incy);
#else
                    dtemp = Cdotu(n, x_ref, incx, y_ref, incy);
#endif
                    ctemp = Cdotu(n, x, incx, y, incy);
                    REAL_REF diff = abs(ctemp - dtemp);
                    if (diff > EPSILON) {
#if defined VERBOSE_TEST
                        printf("error: ");
                        printnum(diff);
                        printf("\n");
#endif
                        errorflag = TRUE;
                    }
                    if (maxdiff < diff)
                        maxdiff = diff;
                    j++;
                }
                delete[] x;
                delete[] y;
                delete[] x_ref;
                delete[] y_ref;
            }
        }
    }
    if (errorflag == TRUE) {
        printf("error: ");
        printnum(maxdiff);
        printf("\n");
        printf("*** Testing Cdotu failed ***\n");
        exit(1);
    } else {
        printf("maxerror: ");
        printnum(maxdiff);
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Cdotu start ***\n");
    Cdotu_test();
    printf("*** Testing Cdotu successful ***\n");
    return (0);
}
