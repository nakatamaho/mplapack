/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Claset.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_LDA 15
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Claset_test2(const char *uplo) {
    int errorflag = FALSE;
    int j = 0;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int m = MIN_M; m < MAX_M; m++) {
            for (int lda = max(1, m); lda < MAX_LDA; lda++) {
                COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
                COMPLEX_REF alpha_ref, beta_ref;
                COMPLEX *A = new COMPLEX[matlen(lda, n)];
                COMPLEX alpha, beta;

                j = 0;
#if defined VERBOSE_TEST
                printf("#n:%d m:%d lda: %d uplo %s\n", n, m, lda, uplo);
#endif
                while (j < MAX_ITER) {
                    set_random_vector(A_ref, A, matlen(lda, n));
                    set_random_number(alpha_ref, alpha);
                    set_random_number(beta_ref, beta);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    zlaset_f77(uplo, &m, &n, &alpha_ref, &beta_ref, A_ref, &lda);
#else
                    Claset(uplo, m, n, alpha_ref, beta_ref, A_ref, lda);
#endif
                    Claset(uplo, m, n, alpha, beta, A, lda);

                    diff = infnorm(A_ref, A, matlen(lda, n), 1);
                    if (diff > EPSILON) {
                        errorflag = TRUE;
                        printf("Error\n");
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
                delete[] A_ref;
            }
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Claset failed ***\n");
        exit(1);
    }
}

void Claset_test(void) {
    Claset_test2("U");
    Claset_test2("L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Claset start ***\n");
    Claset_test();
    printf("*** Testing Claset successful ***\n");
    return (0);
}
