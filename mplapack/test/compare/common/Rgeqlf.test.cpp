/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rsytrd.debug.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 1
#define MAX_N 30
#define MIN_M 1
#define MAX_M 30
#define MAX_LDA 30
#define MAX_ITER 2

REAL_REF maxdiff = 0.0;

void Rgeqlf_test() {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    INTEGER info, lwork;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int m = MIN_M; m < MAX_M; m++) {
            for (int lda = max(m, 1); lda < MAX_LDA; lda++) {

                REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
                REAL_REF *tau_ref = new REAL_REF[veclen(min(m, n), 1)];

                REAL *A = new REAL[matlen(lda, n)];
                REAL *tau = new REAL[veclen(min(m, n), 1)];

#if defined VERBOSE_TEST
                printf("n:%d m:%d lda %d\n", n, m, lda);
#endif
                lwork_ref = -1;
                lwork = -1;
                REAL_REF *work_ref = new REAL_REF[1];
                REAL *work = new REAL[1];
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dgeqlf_f77(&m, &n, A_ref, &lda, tau_ref, work_ref, &lwork_ref, &info_ref);
#else
                Rgeqlf(m, n, A_ref, lda, tau_ref, work_ref, lwork_ref, info_ref);
#endif
                Rgeqlf(m, n, A, lda, tau, work, lwork, info);
                lwork_ref = (int)cast2double(work_ref[0]);
                lwork = (int)cast2double(work[0]);
#if defined VERBOSE_TEST
                printf("optimized worksize by Rgeqlf %d : by dgeqlf %d.\n", (int)lwork, (int)lwork_ref);
#endif
                delete[] work;
                delete[] work_ref;
                lwork_ref = max(lwork_ref, (INTEGER_REF)1);
                lwork = max(lwork, (INTEGER)1);
                work_ref = new REAL_REF[lwork_ref];
                work = new REAL[lwork];
                j = 0;
                while (j < MAX_ITER) {
                    set_random_vector(A_ref, A, matlen(lda, n));
                    set_random_vector(tau_ref, tau, veclen(min(m, n), 1));
                    set_random_vector(work_ref, work, veclen(lwork, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                    dgeqlf_f77(&m, &n, A_ref, &lda, tau_ref, work_ref, &lwork_ref, &info_ref);
#else
                    Rgeqlf(m, n, A_ref, lda, tau_ref, work_ref, lwork_ref, info_ref);
#endif
                    Rgeqlf(m, n, A, lda, tau, work, lwork, info);

                    if (info != info_ref) {
                        printf("info differ! %d, %d\n", (int)info, (int)info_ref);
                        errorflag = TRUE;
                    }
                    if (info < 0) {
                        j++;
                        continue;
                    }
                    diff = infnorm(A_ref, A, matlen(lda, n), 1);
                    if (diff > EPSILON2) {
                        printf("A = ");
                        printmat(n, n, A, lda);
                        printf("\n");
                        printf("error in A: ");
                        printnum(diff);
                        printf("\n");
                        errorflag = TRUE;
                    }
                    if (maxdiff < diff)
                        maxdiff = diff;

                    diff = infnorm(tau_ref, tau, veclen(min(n, m), 1), 1);
                    if (diff > EPSILON2) {
                        printf("error in tau: ");
                        printnum(diff);
                        printf("\n");
                        errorflag = TRUE;
                    }
                    if (maxdiff < diff)
                        maxdiff = diff;

                    diff = infnorm(work_ref, work, veclen(lwork, 1), 1);
                    if (diff > EPSILON2) {
                        printf("error in work: ");
                        printnum(diff);
                        printf("\n");
                        errorflag = TRUE;
                    }

                    if (maxdiff < diff)
                        maxdiff = diff;
                    j++;
#if defined VERBOSE_TEST
                    printf("max error: ");
                    printnum(maxdiff);
                    printf("\n");
#endif
                }
                delete[] work;
                delete[] work_ref;
                delete[] tau;
                delete[] A;
                delete[] tau_ref;
                delete[] A_ref;
            }
            if (errorflag == TRUE) {
                printf("*** Testing Rgeqlf failed ***\n");
                exit(1);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rgeqlf start ***\n");
    Rgeqlf_test();
    printf("*** Testing Rgeqlf successful ***\n");
    return (0);
}
