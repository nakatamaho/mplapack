/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgesv.debug.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 1
#define MAX_N 10
#define MIN_NRHS 1
#define MAX_NRHS 10
#define MAX_LDA 9
#define MAX_LDB 9
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Rgesv_test() {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;

    for (int n = MIN_N; n <= MAX_N; n++) {
        for (int nrhs = MIN_NRHS; nrhs <= n; nrhs++) {
            for (int lda = max(n, 1); lda <= MAX_LDA; lda++) {
                for (int ldb = max(n, 1); ldb <= MAX_LDB; ldb++) {

                    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
                    REAL_REF *B_ref = new REAL_REF[matlen(ldb, nrhs)];
                    INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(n, 1)];

                    REAL *A = new REAL[matlen(lda, n)];
                    REAL *B = new REAL[matlen(ldb, nrhs)];
                    INTEGER *ipiv = new INTEGER[veclen(n, 1)];
#if defined VERBOSE_TEST
                    printf("#n:%d lda %d nrhs %d ldb %d\n", n, lda, nrhs, ldb);
#endif
                    j = 0;
                    while (j < MAX_ITER) {
                        set_random_vector(A_ref, A, matlen(lda, n));
                        set_random_vector(B_ref, B, matlen(ldb, nrhs));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                        dgesv_f77(&n, &nrhs, A_ref, &lda, ipiv_ref, B_ref, &ldb, &info_ref);
#else
                        Rgesv(n, nrhs, A_ref, lda, ipiv_ref, B_ref, ldb, info_ref);
#endif
                        Rgesv(n, nrhs, A, lda, ipiv, B, ldb, info);

                        if (info < 0) {
                            printf("info %d error\n", -(int)info);
                        }
                        if (info_ref != info) {
                            printf("info differ! %d, %d\n", (int)info_ref, (int)info);
                            errorflag = TRUE;
                        }
                        diff = infnorm(B_ref, B, matlen(ldb, nrhs), 1);
                        if (diff > EPSILON7) {
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
                    delete[] ipiv_ref;
                    delete[] B_ref;
                    delete[] A_ref;
                    delete[] ipiv;
                    delete[] B;
                    delete[] A;
                }
                if (errorflag == TRUE) {
                    printf("*** Testing Rgesv failed ***\n");
                    exit(1);
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rgesv start ***\n");
    Rgesv_test();
    printf("*** Testing Rgesv successful ***\n");
    return (0);
}
