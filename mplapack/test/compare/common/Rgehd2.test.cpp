/*
 * Copyright (c) 2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rpotf2.debug.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N 20   // should not be so large
#define MAX_LDA 20 // should not be so large
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rgehd2_test(void) {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, ilo_ref, ihi_ref;
    INTEGER info, ilo, ihi;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
            REAL_REF *A_ref    = new REAL_REF[matlen(lda, n)];
            REAL_REF *tau_ref  = new REAL_REF[veclen(n - 1, 1)];
            REAL_REF *work_ref = new REAL_REF[veclen(n, 1)];
            REAL *A    = new REAL[matlen(lda, n)];
            REAL *tau  = new REAL[veclen(n - 1, 1)];
            REAL *work = new REAL[veclen(n, 1)];
            ilo = 1;
            ilo_ref = 1;
            ihi = n;
            ihi_ref = n;
#if defined VERBOSE_TEST
            printf("n:%d lda %d \n", n, lda);
#endif
            j = 0;
            while (j < MAX_ITER) {
                set_random_vector(A_ref, A, matlen(lda, n));
#if defined VERBOSE_TEST
                printf("Aorg=");printmat(n,n,A,lda);printf("\n");
                printf("Aorg_ref=");printmat(n,n,A_ref,lda);printf("\n");
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dgehd2_f77(&n, &ilo_ref, &ihi_ref, A_ref, &lda, tau_ref, work_ref, &info_ref);
#else
                Rgehd2(n, ilo_ref, ihi_ref, A_ref, lda, tau_ref, work_ref, info_ref);
#endif
                Rgehd2(n, ilo, ihi, A, lda, tau, work, info);

                if (info < 0) {
                    printf("info %d error\n", -(int)info);
                    errorflag = TRUE;
                }
                if (info_ref != info) {
                    printf("info error! %d, %d\n", (int)info_ref, (int)info);
                    errorflag = TRUE;
                }
#if defined VERBOSE_TEST
                printf("h=");printmat(n,n,A,lda);printf("\n");
                printf("h_ref=");printmat(n,n,A_ref,lda);printf("\n");
                printf("tau=");printvec(tau,n-1);printf("\n");
                printf("tau_ref=");printvec(tau_ref,n-1);printf("\n");
#endif
                diff = infnorm(A_ref, A, matlen(lda, n), 1);
                if (diff > EPSILON2) {
                    printf("n:%d lda %d\n", n, lda);
                    printf("error1: ");
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
                if (errorflag == TRUE) {
                    printf("*** Testing Rgehd2 failed ***\n");
                    exit(1);
                }
            }
            delete[] work;
            delete[] work_ref;
            delete[] tau;
            delete[] tau_ref;
            delete[] A;
            delete[] A_ref;
        }
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rgehd2 start ***\n");
    Rgehd2_test();
    printf("*** Testing Rgehd2 successful ***\n");
    return (0);
}
