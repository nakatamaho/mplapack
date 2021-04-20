/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlasr.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N 15
#define MIN_M 0
#define MAX_M 15
#define MAX_LDA 17
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rlasr_test2(const char *side, const char *pivot, const char *direct) {
    int errorflag = FALSE;
    int j = 0;
    int cdim = 0;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int m = MIN_M; m < MAX_M; m++) {
            for (int lda = max(m, 1); lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
                printf("# n:%d m:%d lda:%d, side %s, pivot %s direct %s\n", n, m, lda, side, pivot, direct);
#endif
                if (Mlsame(side, "L"))
                    cdim = m - 1;
                if (Mlsame(side, "R"))
                    cdim = n - 1;

                REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
                REAL_REF *c_ref = new REAL_REF[max(cdim, 1)];
                REAL_REF *s_ref = new REAL_REF[max(cdim, 1)];
                REAL *A = new REAL[matlen(lda, n)];
                REAL *c = new REAL[max(cdim, 1)];
                REAL *s = new REAL[max(cdim, 1)];

                j = 0;
                while (j < MAX_ITER) {
                    set_random_vector(A_ref, A, matlen(lda, n));
                    set_random_vector(c_ref, c, max(cdim, 1));
                    set_random_vector(s_ref, s, max(cdim, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    dlasr_f77(side, pivot, direct, &m, &n, c_ref, s_ref, A_ref, &lda);
#else
                    Rlasr(side, pivot, direct, m, n, c_ref, s_ref, A_ref, lda);
#endif
                    Rlasr(side, pivot, direct, m, n, c, s, A, lda);

                    diff = infnorm(A_ref, A, matlen(lda, n), 1);
                    if (diff > EPSILON) {
                        printf("# n:%d m:%d lda:%d, side %s, pivot %s direct %s\n", n, m, lda, side, pivot, direct);
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
                delete[] s_ref;
                delete[] c_ref;
                delete[] A_ref;
                delete[] s;
                delete[] c;
                delete[] A;
            }
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Rlasr failed ***\n");
        exit(1);
    }
}

void Rlasr_test(void) {
    Rlasr_test2("L", "V", "F");
    Rlasr_test2("R", "V", "F");
    Rlasr_test2("L", "T", "F");
    Rlasr_test2("R", "T", "F");
    Rlasr_test2("L", "B", "F");
    Rlasr_test2("R", "B", "F");
    Rlasr_test2("L", "V", "B");
    Rlasr_test2("R", "V", "B");
    Rlasr_test2("L", "T", "B");
    Rlasr_test2("R", "T", "B");
    Rlasr_test2("L", "B", "B");
    Rlasr_test2("R", "B", "B");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rlasr start ***\n");
    Rlasr_test();
    printf("*** Testing Rlasr successful ***\n");
    return (0);
}
