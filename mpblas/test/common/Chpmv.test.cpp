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

#define MIN_INCX -1
#define MAX_INCX 6
#define MIN_INCY -1
#define MAX_INCY 6
#define MIN_N 2
#define MAX_N 10
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

void Chpmv_test2(const char *uplo) {
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
        for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
            for (int n = MIN_N; n < MAX_N; n++) {
#if defined VERBOSE_TEST
                printf("#n is %d, incx is %d, incy is %d, uplo is %s.\n", n, incx, incy, uplo);
#endif
                COMPLEX_REF alpha_ref;
                COMPLEX_REF beta_ref;
                COMPLEX_REF *x_ref;
                COMPLEX_REF *y_ref;
                COMPLEX_REF *AP_ref;

                COMPLEX alpha;
                COMPLEX beta;
                COMPLEX *x;
                COMPLEX *y;
                COMPLEX *AP;

                x_ref = new COMPLEX_REF[veclen(n, incx)];
                y_ref = new COMPLEX_REF[veclen(n, incy)];
                AP_ref = new COMPLEX_REF[vecplen(n)];
                x = new COMPLEX[veclen(n, incx)];
                y = new COMPLEX[veclen(n, incy)];
                AP = new COMPLEX[vecplen(n)];

                for (int i = 0; i < MAX_ITER; i++) {
                    set_random_vector(AP_ref, AP, vecplen(n));
                    set_random_vector(x_ref, x, veclen(n, incx));
                    set_random_vector(y_ref, y, veclen(n, incy));
                    set_random_number(alpha_ref, alpha);
                    set_random_number(beta_ref, beta);

                    mplapack_errno = 0;
                    blas_errno = 0;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    zhpmv_f77(uplo, &n, &alpha_ref, AP_ref, x_ref, &incx, &beta_ref, y_ref, &incy);
                    mplapack_errno1 = blas_errno;
#else
                    Chpmv(uplo, n, alpha_ref, AP_ref, x_ref, incx, beta_ref, y_ref, incy);
                    mplapack_errno1 = mplapack_errno;
#endif
                    Chpmv(uplo, n, alpha, AP, x, incx, beta, y, incy);
                    mplapack_errno2 = mplapack_errno;
#if defined VERBOSE_TEST
                    printf("errno: mplapack %d, ref %d\n", mplapack_errno1, mplapack_errno2);
#endif
                    if (mplapack_errno1 != mplapack_errno2) {
#if defined VERBOSE_TEST
                        printf("error in Mxerbla!!\n");
#endif
                        errorflag = TRUE;
                    }
                    REAL_REF diff = infnorm(y_ref, y, veclen(n, incy), 1);
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
                }
                delete[] AP_ref;
                delete[] y_ref;
                delete[] x_ref;
                delete[] AP;
                delete[] y;
                delete[] x;
            }
        }
    }
    if (errorflag == TRUE) {
        printf("error: ");
        printnum(maxdiff);
        printf("\n");
        printf("*** Testing Chpmv failed ***\n");
        exit(1);
    } else {
        printf("maxerror: ");
        printnum(maxdiff);
        printf("\n");
    }
}

void Chpmv_test() {
    Chpmv_test2("U");
    Chpmv_test2("L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Chpmv start ***\n");
    Chpmv_test();
    printf("*** Testing Chpmv successful ***\n");
    return (0);
}
