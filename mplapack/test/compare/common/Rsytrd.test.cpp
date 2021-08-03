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
#include <mplapack_compare_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N -1
#define MAX_N 30
#define MAX_LDA 30
#define MAX_ITER 2

REAL_REF maxdiff = 0.0;

void Rsytrd_test2(const char *uplo) {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    INTEGER info, lwork;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int lda = max(n, 1); lda < MAX_LDA; lda++) {

            REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *d_ref = new REAL_REF[veclen(n, 1)];
            REAL_REF *e_ref = new REAL_REF[veclen(n - 1, 1)];
            REAL_REF *tau_ref = new REAL_REF[veclen(n - 1, 1)];

            REAL *A = new REAL[matlen(lda, n)];
            REAL *d = new REAL[veclen(n, 1)];
            REAL *e = new REAL[veclen(n - 1, 1)];
            REAL *tau = new REAL[veclen(n - 1, 1)];

#if defined VERBOSE_TEST
            printf("#uplo %s, n:%d lda %d\n", uplo, n, lda);
#endif
            lwork_ref = -1;
            lwork = -1;
            REAL_REF *work_ref = new REAL_REF[1];
            REAL *work = new REAL[1];
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
            dsytrd_f77(uplo, &n, A_ref, &lda, d_ref, e_ref, tau_ref, work_ref, &lwork_ref, &info_ref);
#else
            Rsytrd(uplo, n, A_ref, lda, d_ref, e_ref, tau_ref, work_ref, lwork_ref, info_ref);
#endif
            Rsytrd(uplo, n, A, lda, d, e, tau, work, lwork, info);
            lwork_ref = (int)cast2double(work_ref[0]);
            lwork = (int)cast2double(work[0]);
#if defined VERBOSE_TEST
            printf("optimized worksize by Rsytrd %d : by dsytrd %d.\n", (int)lwork, (int)lwork_ref);
#endif
#ifdef DUMMY
            // comparison of workspace is nonsense...
            if (worksize != worksized)
                printf("error in worksize\n");
#endif
            delete[] work;
            delete[] work_ref;
            work_ref = new REAL_REF[max(1, (int)lwork_ref)];
            work = new REAL[max(1, (int)lwork)];
            j = 0;
            while (j < MAX_ITER) {
                set_random_vector(A_ref, A, matlen(lda, n));
                set_random_vector(d_ref, d, veclen(n, 1));
                set_random_vector(e_ref, e, veclen(n - 1, 1));
                set_random_vector(tau_ref, tau, veclen(n - 1, 1));
                set_random_vector(work_ref, work, veclen(lwork, 1));

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dsytrd_f77(uplo, &n, A_ref, &lda, d_ref, e_ref, tau_ref, work_ref, &lwork_ref, &info_ref);
#else
                Rsytrd(uplo, n, A_ref, lda, d_ref, e_ref, tau_ref, work_ref, lwork_ref, info_ref);
#endif
                Rsytrd(uplo, n, A, lda, d, e, tau, work, lwork, info);

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
                diff = infnorm(d_ref, d, veclen(n, 1), 1);
                if (diff > EPSILON2) {
                    printf("error in d: ");
                    printnum(diff);
                    printf("\n");
                    errorflag = TRUE;
                }
                if (maxdiff < diff)
                    maxdiff = diff;

                diff = infnorm(e_ref, e, veclen(n - 1, 1), 1);
                if (diff > EPSILON2) {
                    printf("error in e: ");
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
            delete[] e;
            delete[] d;
            delete[] A;
            delete[] tau_ref;
            delete[] e_ref;
            delete[] d_ref;
            delete[] A_ref;
        }
        if (errorflag == TRUE) {
            printf("*** Testing Rsytrd failed ***\n");
            exit(1);
        }
    }
}

void Rsytrd_test(void) {
    Rsytrd_test2("U");
    Rsytrd_test2("L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rsytrd start ***\n");
    Rsytrd_test();
    printf("*** Testing Rsytrd successful ***\n");
    return (0);
}
