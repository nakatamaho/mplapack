/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rpocon.debug.cpp,v 1.2 2010/08/19 01:17:55 nakatamaho Exp $
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

void Rpocon_test2(const char *uplo) {
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
            printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
#endif
            REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *work_ref = new REAL_REF[max(1, n * 3)];
            INTEGER_REF *iwork_ref = new INTEGER_REF[max(1, n)];

            REAL *A = new REAL[matlen(lda, n)];
            REAL *work = new REAL[max(1, n * 3)];
            INTEGER *iwork = new INTEGER[max(1, n)];

            j = 0;
            while (j < MAX_ITER) {
                set_random_psdmat(A_ref, A, lda, n);
/* First, calculate norm of matrix*/
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                anorm_ref = dlange_f77("1", &n, &n, A_ref, &lda, work_ref);
#else
                anorm_ref = Rlange("1", n, n, A_ref, lda, work_ref);
#endif
                anorm = Rlange("1", n, n, A, lda, work);
/* second, do Cholesky factorization via Rpotrf */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dpotrf_f77(uplo, &n, A_ref, &lda, &info_ref);
#else
                Rpotrf(uplo, n, A_ref, lda, info_ref);
#endif
                Rpotrf(uplo, n, A, lda, info);
                if (info > 0) {
#if defined VERBOSE_TEST
                    printf("non psd matrix in %d-th (not an error)\n", (int)info);
#endif
                    break;
                }
/* third, calculate condition number */
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dpocon_f77(uplo, &n, A_ref, &lda, &anorm_ref, &rcond_ref, work_ref, iwork_ref, &info_ref);
#else
                Rpocon(uplo, n, A_ref, lda, anorm_ref, rcond_ref, work_ref, iwork_ref, info_ref);
#endif
                Rpocon(uplo, n, A, lda, anorm, rcond, work, iwork, info);
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
                    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
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
            delete[] iwork;
            delete[] work;
            delete[] A;
            delete[] iwork_ref;
            delete[] work_ref;
            delete[] A_ref;
        }
        if (errorflag == TRUE) {
            printf("*** Testing Rpocon failed ***\n");
            exit(1);
        }
    }
}

void Rpocon_test(void) {
    Rpocon_test2("U");
    Rpocon_test2("L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rpocon start ***\n");
    Rpocon_test();
    printf("*** Testing Rpocon successful ***\n");
    return (0);
}
