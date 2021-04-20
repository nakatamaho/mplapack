/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rsyev.debug.cpp,v 1.14 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N 15
#define MAX_LDA 15
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

template <class X> X infnorm(X *vec1, X *vec2, int len, int inc) {
    X ctemp;
    X inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
        ctemp = vec1[i] - vec2[i];
        if (inorm < abs(ctemp)) {
            inorm = abs(ctemp);
        }
    }
    return inorm;
}

template <class X> X check_eigvec(X *A, X *eigmat, X *w, int n, int lda, const char *uplo) {
    X *tmpmat1 = new X[n * n];
    X *tmpmat2 = new X[n * n];
    X *tmpmat3 = new X[n * n];
    X mtemp, diff;

    if (Mlsame(uplo, "U")) {
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                tmpmat1[i + j * n] = A[i + j * lda];
                tmpmat1[j + i * n] = A[i + j * lda];
            }
        }
    }
    if (Mlsame(uplo, "L")) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                tmpmat1[i + j * n] = A[i + j * lda];
                tmpmat1[j + i * n] = A[i + j * lda];
            }
        }
    }

    for (int i = 0; i < n * n; i++) {
        tmpmat2[i] = 0.0;
    }
    for (int i = 0; i < n; i++) {
        tmpmat2[i + i * n] = w[i];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mtemp = 0.0;
            for (int p = 0; p < n; p++) {
                for (int q = 0; q < n; q++) {
                    mtemp = mtemp + eigmat[p + i * lda] * tmpmat1[p + q * n] * eigmat[q + j * lda];
                }
            }
            tmpmat3[i + j * n] = mtemp;
        }
    }
    diff = infnorm(tmpmat2, tmpmat3, n * n, 1);
    delete[] tmpmat3;
    delete[] tmpmat2;
    delete[] tmpmat1;
    return diff;
    /*
      ( U' * A * U ) = W
      printf("A = ");  printmat(n, n, A, n); printf("\n\n");
      printf("W = ");  printmat(n, n, tmpmat2, n); printf("\n\n");
      printf("U = ");  printmat(n, n, eigmat, n); printf("\n\n");
    */
}

void Rsyev_test2(const char *jobz, const char *uplo) {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    INTEGER info, lwork;
    REAL_REF diff = 0.0;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
            REAL_REF *Aorg_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
            REAL_REF *w_ref = new REAL_REF[veclen(n, 1)];
            REAL *Aorg = new REAL[matlen(lda, n)];
            REAL *A = new REAL[matlen(lda, n)];
            REAL *w = new REAL[veclen(n, 1)];
#if defined VERBOSE_TEST
            printf("# jobz %s, uplo %s,  n:%d lda %d\n", jobz, uplo, n, lda);
#endif
            REAL_REF *work_ref = new REAL_REF[1];
            REAL *work = new REAL[1];
            // these workspace query might not be the same value.
            lwork_ref = -1;
            lwork = -1;
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
            dsyev_f77(jobz, uplo, &n, A_ref, &lda, w_ref, work_ref, &lwork_ref, &info_ref);
#else
            Rsyev(jobz, uplo, n, A_ref, lda, w_ref, work_ref, lwork_ref, &info_ref);
#endif
            Rsyev(jobz, uplo, n, A, lda, w, work, lwork, &info);

            lwork_ref = (int)cast2double(work_ref[0]);
            lwork = (int)cast2double(work[0]);
#if defined VERBOSE_TEST
            printf("optimized worksize by dsyev %d : by Rsyev %d.\n", (int)lwork_ref, (int)lwork);
#endif
#ifdef DUMMY
            // comparison of workspace is nonsense...
            if (worksize != worksize_ref)
                printf("error in worksize\n");
#endif
            delete[] work_ref;
            delete[] work;
            work_ref = new REAL_REF[max(1, (int)lwork_ref)];
            work = new REAL[max(1, (int)lwork)];

            j = 0;
            while (j < MAX_ITER) {
                set_random_vector(A_ref, A, matlen(lda, n));
                set_random_vector(w_ref, w, veclen(n, 1));

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                int iOne = 1, _n;
                _n = matlen(lda, n);
                dcopy_f77(&_n, A_ref, &iOne, Aorg_ref, &iOne);
                Rcopy(matlen(lda, n), A, 1, Aorg, 1);
#else
                Rcopy(matlen(lda, n), A_ref, 1, Aorg_ref, 1);
                Rcopy(matlen(lda, n), A, 1, Aorg, 1);
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                dsyev_f77(jobz, uplo, &n, A_ref, &lda, w_ref, work_ref, &lwork_ref, &info_ref);
#else
                Rsyev(jobz, uplo, n, A_ref, lda, w_ref, work_ref, lwork_ref, &info_ref);
#endif
                Rsyev(jobz, uplo, n, A, lda, w, work, lwork, &info);

                if (info < 0) {
                    printf("info %d error\n", -(int)info);
                    errorflag = TRUE;
                    exit(1);
                }
                if (info != info_ref) {
                    printf("info differ! %d, %d\n", (int)info_ref, (int)info);
                    errorflag = TRUE;
                }
                if (Mlsame(jobz, "V")) {
                    diff = check_eigvec(Aorg, A, w, (int)n, (int)lda, uplo);
                }
                if (diff > EPSILON) {
                    printf("error in eigvec1: ");
                    printnum(diff);
                    printf("\n");
                    errorflag = TRUE;
                }
                if (maxdiff < diff)
                    maxdiff = diff;

                if (Mlsame(jobz, "V")) {
                    diff = check_eigvec(Aorg_ref, A_ref, w_ref, (int)n, (int)lda, uplo);
                }
                if (diff > EPSILON) {
                    printf("error in eigvec2: ");
                    printnum(diff);
                    printf("\n");
                    errorflag = TRUE;
                }
                if (maxdiff < diff)
                    maxdiff = diff;

                diff = infnorm(w_ref, w, veclen(n, 1), 1);
                if (diff > EPSILON) {
                    printf("error in w: ");
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
            delete[] w;
            delete[] A;
            delete[] Aorg;
            delete[] w_ref;
            delete[] A_ref;
            delete[] Aorg_ref;
        }
        if (errorflag == TRUE) {
            printf("*** Testing Rsyev failed ***\n");
            exit(1);
        }
    }
}

void Rsyev_test(void) {
    Rsyev_test2("V", "U");
    Rsyev_test2("V", "L");
    Rsyev_test2("N", "U");
    Rsyev_test2("N", "L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Rsyev start ***\n");
    Rsyev_test();
    printf("*** Testing Rsyev successful ***\n");
    return (0);
}
