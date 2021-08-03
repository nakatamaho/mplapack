/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlatrd.debug.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 4
#define MAX_N 10
#define MIN_LDA 4
#define MAX_LDA 10
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Rlatrs_test2(const char *uplo, const char *trans, const char *diag, const char *normin) {
    int errorflag = FALSE;
    REAL_REF scale_ref;
    REAL scale;
    INTEGER_REF info_ref;
    INTEGER info;
    INTEGER iter;
    REAL_REF diff;
    INTEGER_REF n, lda;
    REAL rtmp;
    REAL_REF rtmp_ref;

    for (n = MIN_N; n <= MAX_N; n++) {
        for (lda = max(1, (int)n); lda <= MAX_LDA; lda++) {
#if defined VERBOSE_TEST
            printf("# n %d, lda %d, uplo %s, trans %s, diag %s, normin %s\n", (int)n, (int)lda, uplo, trans, diag, normin);
#endif
            REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *B_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *x_ref = new REAL_REF[veclen(n, 1)];
            REAL_REF *cnorm_ref = new REAL_REF[veclen(n, 1)];
            REAL *A = new REAL[matlen(lda, n)];
            REAL *x = new REAL[veclen(n, 1)];
            REAL *cnorm = new REAL[veclen(n, 1)];

            REAL_REF *y_ref = new REAL_REF[veclen(n, 1)];
            REAL *y = new REAL[veclen(n, 1)];

            for (iter = 0; iter < MAX_ITER; iter++) {
                set_random_vector(A_ref, A, matlen(lda, n));
                set_random_vector(A_ref, A, matlen(lda, n));
                set_random_vector(x_ref, x, veclen(n, 1));
                set_random_vector(cnorm_ref, cnorm, veclen(n, 1));
                /* A*x=b; b would be very large. so we choose x as the answer. */
                for (int p = 0; p < n; p++) {
                    for (int q = 0; q < n; q++) {
                        B_ref[p + q * lda] = 0.0;
                    }
                }
                if (Mlsame(uplo, "U")) {
                    for (int p = 0; p < n; p++) {
                        for (int q = p; q < n; q++) {
                            B_ref[p + q * lda] = A_ref[p + q * lda];
                        }
                    }
                }
                if (Mlsame(uplo, "L")) {
                    for (int p = 0; p < n; p++) {
                        for (int q = 0; q <= p; q++) {
                            B_ref[p + q * lda] = A_ref[p + q * lda];
                        }
                    }
                }
                if (Mlsame(diag, "U")) {
                    for (int p = 0; p < n; p++) {
                        B_ref[p + p * lda] = 1.0;
                    }
                }
                // printf("A"); printmat(n, n, B_ref, lda); printf("\n");
                for (int p = 0; p < n; p++) {
                    rtmp = 0.0;
                    rtmp_ref = 0.0;
                    for (int q = 0; q < n; q++) {
                        if (Mlsame(trans, "N")) {
                            rtmp_ref += B_ref[p + q * lda] * x_ref[q];
                        } else {
                            rtmp_ref += B_ref[q + p * lda] * x_ref[q];
                        }
                    }
                    y_ref[p] = rtmp_ref;
                }
                // printf("ans"); printvec(x_ref, veclen(n, 1)); printf("\n");
                // printf("A"); printmat(n, n, B_ref, lda); printf("\n");
                for (int p = 0; p < n; p++) {
                    x_ref[p] = y_ref[p];
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    x[p] = y_ref[p];
#elif defined ___MPLAPACK_BUILD_WITH_GMP___
                    x[p] = cast2mpf_class(y_ref[p]);
#elif defined ___MPLAPACK_BUILD_WITH_QD___
                    x[p] = cast2qd_real(y_ref[p]);
#elif defined ___MPLAPACK_BUILD_WITH_DD___
                    x[p] = cast2dd_real(y_ref[p]);
#elif defined ___MPLAPACK_BUILD_WITH_DOUBLE___
                    x[p] = cast2double(y_ref[p]);
#elif defined ___MPLAPACK_BUILD_WITH__FLOAT128___
                    x[p] = cast2_Float128(y_ref[p]);
#endif
                }
                // printf("y_ref"); printvec(y_ref, veclen(n, 1)); printf("\n");

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dlatrs_f77(uplo, trans, diag, normin, &n, A_ref, &lda, x_ref, &scale_ref, cnorm_ref, &info_ref);
#else
                Rlatrs(uplo, trans, diag, normin, n, A_ref, lda, x_ref, scale_ref, cnorm_ref, info_ref);
#endif
                Rlatrs(uplo, trans, diag, normin, n, A, lda, x, scale, cnorm, info);
                if (info != info_ref) {
                    printf("info differ! %d, %d\n", (int)info, (int)info_ref);
                    errorflag = TRUE;
                }
                diff = infnorm(x_ref, x, veclen(n, 1), 1);
                if (diff > EPSILON10) {
                    printf("error in x: ");
                    printnum(diff);
                    printf("\n");
                    printf("x_ref");
                    printvec(x_ref, veclen(n, 1));
                    printf("\n");
                    printf("x");
                    printvec(x, veclen(n, 1));
                    printf("\n");
                    errorflag = TRUE;
                }
                if (maxdiff < diff)
                    maxdiff = diff;

                diff = infnorm(cnorm_ref, cnorm, veclen(n, 1), 1);
                if (diff > EPSILON10) {
                    printf("error in cnorm: ");
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
            delete[] cnorm;
            delete[] x;
            delete[] A;
            delete[] cnorm_ref;
            delete[] x_ref;
            delete[] B_ref;
            delete[] A_ref;
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rlatrs failed ***\n");
        exit(1);
    }
}

void Rlatrs_test() {
    Rlatrs_test2("U", "N", "N", "Y");
    Rlatrs_test2("U", "N", "N", "N");
    Rlatrs_test2("U", "N", "U", "Y");
    Rlatrs_test2("U", "N", "U", "N");

    Rlatrs_test2("U", "T", "N", "Y");
    Rlatrs_test2("U", "T", "N", "N");
    Rlatrs_test2("U", "T", "U", "Y");
    Rlatrs_test2("U", "T", "U", "N");

    Rlatrs_test2("U", "C", "N", "Y");
    Rlatrs_test2("U", "C", "N", "N");
    Rlatrs_test2("U", "C", "U", "Y");
    Rlatrs_test2("U", "C", "U", "N");

    Rlatrs_test2("L", "N", "N", "Y");
    Rlatrs_test2("L", "N", "N", "N");
    Rlatrs_test2("L", "N", "U", "Y");
    Rlatrs_test2("L", "N", "U", "N");

    Rlatrs_test2("L", "T", "N", "Y");
    Rlatrs_test2("L", "T", "N", "N");
    Rlatrs_test2("L", "T", "U", "Y");
    Rlatrs_test2("L", "T", "U", "N");

    Rlatrs_test2("L", "C", "N", "Y");
    Rlatrs_test2("L", "C", "N", "N");
    Rlatrs_test2("L", "C", "U", "Y");
    Rlatrs_test2("L", "C", "U", "N");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rlatrs start ***\n");
    Rlatrs_test();
    printf("*** Testing Rlatrs successful ***\n");
    return (0);
}
