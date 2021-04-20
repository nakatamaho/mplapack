/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cung2l.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MIN_M 0
#define MAX_M 10
#define MIN_K 0
#define MAX_K 10
#define MIN_LDA 4
#define MAX_LDA 10
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Cung2l_test() {
    int errorflag = FALSE;
    int iter;
    int m, n, k, lda;
    INTEGER_REF info_ref;
    REAL_REF diff;
    INTEGER info;

    for (m = MIN_M; m <= MAX_M; m++) {
        for (n = MIN_N; n <= m; n++) {
            for (k = MIN_K; k <= n; k++) {
                for (lda = max(1, m); lda <= MAX_LDA; lda++) {
#if defined VERBOSE_TEST
                    printf("# m %d n %d k %d lda %d\n", (int)m, (int)n, (int)k, (int)lda);
#endif
                    COMPLEX *A = new COMPLEX[matlen(lda, n)];
                    COMPLEX *tau = new COMPLEX[veclen(k, 1)];
                    COMPLEX *work = new COMPLEX[veclen(n, 1)];

                    COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
                    COMPLEX_REF *tau_ref = new COMPLEX_REF[veclen(k, 1)];
                    COMPLEX_REF *work_ref = new COMPLEX_REF[veclen(n, 1)];

                    for (iter = 0; iter < MAX_ITER; iter++) {
                        set_random_vector(A_ref, A, matlen(lda, n));
                        set_random_vector(tau_ref, tau, veclen(k, 1));
                        set_random_vector(work_ref, work, veclen(n, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                        zung2l_f77(&m, &n, &k, A_ref, &lda, tau_ref, work_ref, &info_ref);
#else
                        Cung2l(m, n, k, A_ref, lda, tau_ref, work_ref, info_ref);
#endif
                        Cung2l(m, n, k, A, lda, tau, work, info);

                        diff = infnorm(A_ref, A, matlen(lda, n), 1);
                        if (diff > EPSILON2) {
                            printf("error in A: ");
                            printnum(diff);
                            printf("\n");
                            errorflag = TRUE;
                        }
                        if (maxdiff < diff)
                            maxdiff = diff;
                        diff = infnorm(tau_ref, tau, veclen(k, 1), 1);
                        if (diff > EPSILON2) {
                            printf("error in t: ");
                            printnum(diff);
                            printf("\n");
                            errorflag = TRUE;
                        }
                        if (maxdiff < diff)
                            maxdiff = diff;

                        diff = infnorm(work_ref, work, veclen(n, 1), 1);
                        if (diff > EPSILON2) {
                            printf("error in t: ");
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
                    delete[] tau_ref;
                    delete[] work_ref;
                    delete[] A_ref;
                    delete[] tau;
                    delete[] work;
                    delete[] A;
                }
            }
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Cung2l failed ***\n");
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Cung2l start ***\n");
    Cung2l_test();
    printf("*** Testing Cung2l successful ***\n");
    return (0);
}
