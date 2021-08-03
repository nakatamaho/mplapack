/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgeequ.debug.cpp,v 1.2 2010/08/19 01:17:55 nakatamaho Exp $
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

#define MIN_N 2
#define MIN_M 2
#define MAX_N 10
#define MAX_M 10
#define MAX_LDA 10
#define MAX_ITER 10
REAL_REF maxdiff = 0.0;

void Rgeequ_test(void) {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;
    REAL_REF rtmp_ref;
    REAL rtmp;

    REAL_REF rowcnd_ref, colcnd_ref, amax_ref;
    REAL rowcnd, colcnd, amax;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int m = MIN_M; m < MAX_M; m++) {
            for (int lda = max(m, 1); lda < MAX_LDA; lda++) {
                REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
                REAL_REF *R_ref = new REAL_REF[veclen(m, 1)];
                REAL_REF *C_ref = new REAL_REF[veclen(n, 1)];

                REAL *A = new REAL[matlen(lda, n)];
                REAL *R = new REAL[veclen(m, 1)];
                REAL *C = new REAL[veclen(n, 1)];

#if defined VERBOSE_TEST
                printf("n:%d m:%d lda %d\n", n, m, lda);
#endif
                j = 0;
                while (j < MAX_ITER) {
                    set_random_vector(A_ref, A, matlen(lda, n));
                    set_random_vector(R_ref, R, veclen(m, 1));
                    set_random_vector(C_ref, C, veclen(n, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    dgeequ_f77(&m, &n, A_ref, &lda, R_ref, C_ref, &rowcnd_ref, &colcnd_ref, &amax_ref, &info_ref);
#else
                    Rgeequ(m, n, A_ref, lda, R_ref, C_ref, rowcnd_ref, colcnd_ref, amax_ref, info_ref);
#endif
                    Rgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info);

                    if (info != info_ref) {
                        printf("Diff info: n:%d m:%d lda %d\n", n, m, lda);
                        errorflag = TRUE;
                    }

                    diff = infnorm(R_ref, R, veclen(m, 1), 1);
                    if (diff > EPSILON) {
                        printf("Diff R: n:%d m:%d lda %d\n", n, m, lda);
                        printnum(diff);
                        printf("\n");
                        printf("R_ref = ");
                        printvec(R_ref, veclen(m, 1));
                        printf("\n");
                        printf("R = ");
                        printvec(R, veclen(m, 1));
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
                    diff = infnorm(C_ref, C, veclen(n, 1), 1);
                    if (diff > EPSILON) {
                        printf("Diff C: n:%d m:%d lda %d\n", n, m, lda);
                        errorflag = TRUE;
                    }
                    if (maxdiff < diff)
                        maxdiff = diff;
#if defined VERBOSE_TEST
                    printf("max error: ");
                    printnum(maxdiff);
                    printf("\n");
                    printf("rowcnd_ref");
                    printnum(rowcnd_ref);
                    printf("\n");
#endif
                    diff = abs(rowcnd_ref - rowcnd);
                    if (diff > EPSILON) {
                        printf("Diff rowcnd: n:%d m:%d lda %d\n", n, m, lda);
                        errorflag = TRUE;
                    }
                    if (maxdiff < diff)
                        maxdiff = diff;
#if defined VERBOSE_TEST
                    printf("max error: ");
                    printnum(maxdiff);
                    printf("\n");
                    printf("colcnd_ref");
                    printnum(colcnd_ref);
                    printf("\n");
#endif
                    diff = abs(colcnd_ref - colcnd);
                    if (diff > EPSILON) {
                        printf("Diff colcnd: n:%d m:%d lda %d\n", n, m, lda);
                        errorflag = TRUE;
                    }
                    if (maxdiff < diff)
                        maxdiff = diff;
#if defined VERBOSE_TEST
                    printf("max error: ");
                    printnum(maxdiff);
                    printf("\n");
#endif
                    diff = abs(amax_ref - amax);
                    if (diff > EPSILON) {
                        printf("Diff amax: n:%d m:%d lda %d\n", n, m, lda);
                        printf("amax_ref = ");
                        printnum(amax_ref);
                        printf("\n");
                        printf("amax = ");
                        printnum(amax);
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
                delete[] A;
                delete[] R;
                delete[] C;
                delete[] A_ref;
                delete[] R_ref;
                delete[] C_ref;
            }
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rgeequ failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rgeequ start ***\n");
    Rgeequ_test();
    printf("*** Testing Rgeequ successful ***\n");
    return (0);
}
