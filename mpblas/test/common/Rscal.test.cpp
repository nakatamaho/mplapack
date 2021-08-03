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

#define MIN_INCX -10
#define MAX_INCX 10
#define MAX_N 100
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Rscal_test() {
    int errorflag = FALSE;
    int mplapack_errno1, mplapack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
        for (int n = 0; n < MAX_N; n++) {
#if defined VERBOSE_TEST
            printf("# n:%d incx:%d\n", n, incx);
#endif
            REAL_REF *x_ref = new REAL_REF[veclen(n, incx)];
            REAL *x = new REAL[veclen(n, incx)];
            int j = 0;
            while (j < MAX_ITER) {
                REAL_REF alpha_ref;
                REAL alpha;

                set_random_number(alpha_ref, alpha);
                set_random_vector(x_ref, x, veclen(n, incx));

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dscal_f77(&n, &alpha_ref, x_ref, &incx);
                mplapack_errno1 = blas_errno;
#else
                Rscal(n, alpha_ref, x_ref, incx);
                mplapack_errno1 = mplapack_errno;
#endif
                Rscal(n, alpha, x, incx);
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
                REAL_REF diff = infnorm(x_ref, x, veclen(n, incx), 1);
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
            delete[] x_ref;
            delete[] x;
        }
    }
    if (errorflag == TRUE) {
        printf("error: ");
        printnum(maxdiff);
        printf("\n");
        printf("*** Testing Rscal failed ***\n");
        exit(1);
    } else {
        printf("maxerror: ");
        printnum(maxdiff);
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rscal start ***\n");
    Rscal_test();
    printf("*** Testing Rscal successful ***\n");
    return (0);
}
