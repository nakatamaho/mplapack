/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgecon.debug.cpp,v 1.2 2010/08/19 01:17:55 nakatamaho Exp $
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

#define MIN_N 3
#define MAX_N 20
#define MAX_LDA 20
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

void Rgecon_test2(const char *norm) {
    int errorflag = FALSE;
    int j = 0;
    REAL_REF anorm_ref;
    REAL_REF rcond_ref;
    REAL_REF diff;
    INTEGER_REF info_ref;
    REAL anorm;
    REAL rcond;
    INTEGER info;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
            printf("n:%d lda %d, norm %s\n", n, lda, norm);
#endif
            REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *work_ref = new REAL_REF[max(1, n * 4)];
            INTEGER_REF *iwork_ref = new INTEGER_REF[max(1, n)];
            INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(n, 1)];

            REAL *A = new REAL[matlen(lda, n)];
            REAL *work = new REAL[max(1, n * 4)];
            INTEGER *iwork = new INTEGER[max(1, n)];
            INTEGER *ipiv = new INTEGER[veclen(n, 1)];

            j = 0;
            while (j < MAX_ITER) {
                set_random_vector(A_ref, A, matlen(lda, n));
/* First, calculate norm of matrix*/
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                anorm_ref = dlange_f77(norm, &n, &n, A_ref, &lda, work_ref);
#else
                anorm_ref = Rlange(norm, n, n, A_ref, lda, work_ref);
#endif
                anorm = Rlange(norm, n, n, A, lda, work);
/* second, do LU factorization vir Rgetrf */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dgetrf_f77(&n, &n, A_ref, &lda, ipiv_ref, &info_ref);
#else
                Rgetrf(n, n, A_ref, lda, ipiv_ref, &info_ref);
#endif
                Rgetrf(n, n, A, lda, ipiv, &info);
/* third, calculate condition number */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dgecon_f77(norm, &n, A_ref, &lda, &anorm_ref, &rcond_ref, work_ref, iwork_ref, &info_ref);
#else
                Rgecon(norm, n, A_ref, lda, anorm_ref, &rcond_ref, work_ref, iwork_ref, &info_ref);
#endif
                Rgecon(norm, n, A, lda, anorm, &rcond, work, iwork, &info);
                if (info_ref != info) {
                    printf("info differ! %d, %d\n", (int)info_ref, (int)info);
                    errorflag = TRUE;
                }
                diff = (rcond_ref - rcond);
#if defined VERBOSE_TEST
                printf("reciprocal to cond num:");
                printnum(rcond_ref);
                printf("\n");
                printf("reciprocal to cond num:");
                printnum(rcond);
                printf("\n");
#endif
                if (diff > EPSILON) {
                    printf("n:%d lda %d, norm %s\n", n, lda, norm);
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
            delete[] ipiv;
            delete[] iwork;
            delete[] work;
            delete[] A;
            delete[] ipiv_ref;
            delete[] iwork_ref;
            delete[] work_ref;
            delete[] A_ref;
        }
        if (errorflag == TRUE) {
            printf("*** Testing Rgecon failed ***\n");
            exit(1);
        }
    }
}

void Rgecon_test(void) {
    Rgecon_test2("1");
    Rgecon_test2("I");
    Rgecon_test2("O");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rgecon start ***\n");
    Rgecon_test();
    printf("*** Testing Rgecon successful ***\n");
    return (0);
}
