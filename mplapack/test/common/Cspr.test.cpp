/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cspr.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
#include <blas.h>
#include <lapack.h>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_INCX -2
#define MAX_INCX 3
#define MIN_N 2
#define MAX_N 40
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Cspr_test2(const char *uplo) {
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
        for (int n = MIN_N; n <= MAX_N; n++) {
#if defined VERBOSE_TEST
            printf("#n is %d, incx is %d ", n, incx);
            printf("uplo is %s \n", uplo);
#endif
            COMPLEX_REF *x_ref;
            COMPLEX_REF *AP_ref;
            COMPLEX_REF alpha_ref;

            COMPLEX *x;
            COMPLEX *AP;
            COMPLEX alpha;

            x_ref = new COMPLEX_REF[veclen(n, incx)];
            AP_ref = new COMPLEX_REF[vecplen(n)];
            x = new COMPLEX[veclen(n, incx)];
            AP = new COMPLEX[vecplen(n)];

            for (int i = 0; i < MAX_ITER; i++) {
                set_random_vector(AP_ref, AP, vecplen(n));
                set_random_vector(x_ref, x, veclen(n, incx));
                set_random_number(alpha_ref, alpha);

                mplapack_errno = 0;
                blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                zspr_f77(uplo, &n, &alpha_ref, x_ref, &incx, AP_ref);
                mplapack_errno1 = blas_errno;
#else
                Cspr(uplo, n, alpha_ref, x_ref, incx, AP_ref);
                mplapack_errno1 = mplapack_errno;
#endif
                Cspr(uplo, n, alpha, x, incx, AP);
                mplapack_errno2 = mplapack_errno;
#if defined VERBOSE_TEST
                printf("errno: mplapack %d, ref %d\n", mplapack_errno1, mplapack_errno2);
#endif
                if (mplapack_errno1 != mplapack_errno2) {
                    printf("error in Mxerbla!!\n");
                }
                REAL_REF diff = infnorm(AP_ref, AP, vecplen(n), 1);
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
            }
            delete[] AP_ref;
            delete[] x_ref;
            delete[] AP;
            delete[] x;
        }
    }
    if (errorflag) {
        printf("*** Testing Cspr failed ***\n");
        exit(1);
    }
}

void Cspr_test() {
    Cspr_test2("U");
    Cspr_test2("L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Cspr start ***\n");
    Cspr_test();
    printf("*** Testing Cspr successful ***\n");
    return (0);
}
