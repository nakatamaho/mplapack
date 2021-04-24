/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cgetri.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N 20
#define MAX_LDA 20
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

void Cgetri_test() {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    REAL_REF diff;
    INTEGER info, lwork;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
            COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
            INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(n, 1)];

            COMPLEX *A = new COMPLEX[matlen(lda, n)];
            INTEGER *ipiv = new INTEGER[veclen(n, 1)];
#if defined VERBOSE_TEST
            printf("#n:%d lda %d\n", n, lda);
#endif
            // these workspace query might not be the same value.
            lwork_ref = -1;
            lwork = -1;
            COMPLEX_REF *work_ref = new COMPLEX_REF[1];
            COMPLEX *work = new COMPLEX[1];
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
            zgetri_f77(&n, A_ref, &lda, ipiv_ref, work_ref, &lwork_ref, &info_ref);
#else
            Cgetri(n, A_ref, lda, ipiv_ref, work_ref, lwork_ref, info_ref);
#endif
            Cgetri(n, A, lda, ipiv, work, lwork, info);

            lwork_ref = (int)cast2double(work_ref[0].real());
            lwork = (int)cast2double(work[0].real());
#if defined VERBOSE_TEST
            printf("optimized worksize by Cgetri %d : by dgetri %d.\n", (int)lwork_ref, (int)lwork);
#endif
#ifdef DUMMY
            // comparison of workspace is nonsense...
            if (worksize != worksize_ref)
                printf("error in worksize\n");
#endif
            delete[] work;
            delete[] work_ref;
            work_ref = new COMPLEX_REF[max(1, (int)lwork_ref)];
            work = new COMPLEX[max(1, (int)lwork)];
            j = 0;
            while (j < MAX_ITER) {
                set_random_vector(A_ref, A, matlen(lda, n));
                set_random_vector(work_ref, work, veclen(lwork, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                zgetrf_f77(&n, &n, A_ref, &lda, ipiv_ref, &info_ref);
                zgetri_f77(&n, A_ref, &lda, ipiv_ref, work_ref, &lwork_ref, &info_ref);
#else
                Cgetrf(n, n, A_ref, lda, ipiv_ref, info_ref);
                Cgetri(n, A_ref, lda, ipiv_ref, work_ref, lwork_ref, info_ref);
#endif
                Cgetrf(n, n, A, lda, ipiv, info);
                Cgetri(n, A, lda, ipiv, work, lwork, info);

                if (info < 0) {
                    printf("info %d error\n", -(int)info);
                }
                if (info_ref != info) {
                    printf("info differ! %d, %d\n", (int)info_ref, (int)info);
                    errorflag = TRUE;
                }
                diff = infnorm(A_ref, A, matlen(lda, n), 1);
                if (diff > EPSILON2) {
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
            delete[] work;
            delete[] work_ref;
            delete[] ipiv_ref;
            delete[] A_ref;
            delete[] ipiv;
            delete[] A;
        }
        if (errorflag == TRUE) {
            printf("*** Testing Cgetri failed ***\n");
            exit(1);
        }
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing Cgetri start ***\n");
    Cgetri_test();
    printf("*** Testing Cgetri successful ***\n");
    return (0);
}
