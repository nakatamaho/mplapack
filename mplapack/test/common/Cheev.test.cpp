/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cheev.debug.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 2
#define MAX_N 10
#define MAX_LDA 10
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

template <class X, class Y> X infnorm(Y *vec1, Y *vec2, int len, int inc) {
    X inorm = 0.0;
    Y ctemp;

    for (int i = 0; i < len * inc; i = i + inc) {
        ctemp = vec1[i] - vec2[i];
        if (inorm < abs(ctemp)) {
            inorm = abs(ctemp);
        }
    }
    return inorm;
}

template <class X, class Y> X check_eigvec(Y *A, Y *eigmat, X *w, int n, int lda, const char *uplo) {
    X diff;
    Y *tmpmat1 = new Y[n * n];
    Y *tmpmat2 = new Y[n * n];
    Y *tmpmat3 = new Y[n * n];
    Y mtemp;

    if (Mlsame(uplo, "U")) {
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                tmpmat1[i + j * n] = A[i + j * lda];
                tmpmat1[j + i * n] = conj(A[i + j * lda]);
            }
        }
        for (int i = 0; i < n; i++) {
            tmpmat1[i + i * n] = A[i + i * lda].real();
        }
    }
    if (Mlsame(uplo, "L")) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                tmpmat1[i + j * n] = A[i + j * lda];
                tmpmat1[j + i * n] = conj(A[i + j * lda]);
            }
        }
        for (int i = 0; i < n; i++) {
            tmpmat1[i + i * n] = A[i + i * lda].real();
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
                    mtemp = mtemp + conj(eigmat[p + i * lda]) * tmpmat1[p + q * n] * eigmat[q + j * lda];
                }
            }
            tmpmat3[i + j * n] = mtemp;
        }
    }
    diff = infnorm<X, Y>(tmpmat2, tmpmat3, n * n, 1);
    delete[] tmpmat3;
    delete[] tmpmat2;
    delete[] tmpmat1;
    return diff;
    /*
      ( U' * A * U ) = W
        printf("A = "); printmat(n, n, tmpmat1, n); printf("\n\n");
        printf("W = "); printmat(n, n, tmpmat2, n); printf("\n\n");
        printf("U = "); printmat(n, n, eigmat, n); printf("\n\n");
    */
}

void Cheev_test2(const char *jobz, const char *uplo) {
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    INTEGER info, lwork;
    REAL_REF diff = 0.0;

    for (int n = MIN_N; n < MAX_N; n++) {
        for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
            COMPLEX_REF *Aorg_ref = new COMPLEX_REF[matlen(lda, n)];
            COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
            REAL_REF *w_ref = new REAL_REF[veclen(n, 1)];
            REAL_REF *rwork_ref = new REAL_REF[max(1, 3 * n - 2)];

            COMPLEX *Aorg = new COMPLEX[matlen(lda, n)];
            COMPLEX *A = new COMPLEX[matlen(lda, n)];
            REAL *w = new REAL[veclen(n, 1)];
            REAL *rwork = new REAL[max(1, 3 * n - 2)];
#if defined VERBOSE_TEST
            printf("# jobz %s, uplo %s,  n:%d lda %d\n", jobz, uplo, n, lda);
#endif
            COMPLEX_REF *work_ref = new COMPLEX_REF[1];
            COMPLEX *work = new COMPLEX[1];
            // these workspace query might not be the same value.
            lwork_ref = -1;
            lwork = -1;

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
            zheev_f77(jobz, uplo, &n, A_ref, &lda, w_ref, work_ref, &lwork_ref, rwork_ref, &info_ref);
#else
            Cheev(jobz, uplo, n, A_ref, lda, w_ref, work_ref, lwork_ref, rwork_ref, info_ref);
#endif
            Cheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info);

            lwork_ref = (int)cast2double(work_ref[0].real());
            lwork = (int)cast2double(work[0].real());
#if defined VERBOSE_TEST
            printf("optimized worksize by zheev %d : by Cheev %d.\n", (int)lwork_ref, (int)lwork);
#endif
#ifdef DUMMY
            // comparison of workspace is nonsense...
            if (worksize != worksize_ref)
                printf("error in worksize\n");
#endif
            delete[] work_ref;
            delete[] work;

            work_ref = new COMPLEX_REF[max(1, (int)lwork_ref)];
            work = new COMPLEX[max(1, (int)lwork)];
            j = 0;
            while (j < MAX_ITER) {
                set_random_vector(A_ref, A, matlen(lda, n));
                set_random_vector(w_ref, w, veclen(n, 1));
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                int iOne = 1, _n;
                _n = matlen(lda, n);
                zcopy_f77(&_n, A_ref, &iOne, Aorg_ref, &iOne);
                Ccopy(matlen(lda, n), A, 1, Aorg, 1);
#else
                Ccopy(matlen(lda, n), A_ref, 1, Aorg_ref, 1);
                Ccopy(matlen(lda, n), A, 1, Aorg, 1);
#endif
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                zheev_f77(jobz, uplo, &n, A_ref, &lda, w_ref, work_ref, &lwork_ref, rwork_ref, &info_ref);
#else
                Cheev(jobz, uplo, n, A_ref, lda, w_ref, work_ref, lwork_ref, rwork_ref, info_ref);
#endif
                Cheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info);
                if (info < 0) {
                    printf("info %d error\n", -(int)info);
                    errorflag = TRUE;
                }
                if (info != info_ref) {
                    printf("info differ! %d, %d\n", (int)info_ref, (int)info);
                    errorflag = TRUE;
                }
                if (Mlsame(jobz, "V")) {
                    diff = check_eigvec<REAL, COMPLEX>(Aorg, A, w, (int)n, (int)lda, uplo);
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
                    diff = check_eigvec<REAL_REF, COMPLEX_REF>(Aorg_ref, A_ref, w_ref, (int)n, (int)lda, uplo);
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
#if defined VERBOSE_TEST
                printf("max error: ");
                printnum(maxdiff);
                printf("\n");
#endif
                j++;
            }
            delete[] work;
            delete[] work_ref;
            delete[] rwork;
            delete[] w;
            delete[] A;
            delete[] Aorg;
            delete[] rwork_ref;
            delete[] w_ref;
            delete[] A_ref;
            delete[] Aorg_ref;
        }
        if (errorflag == TRUE) {
            printf("*** Testing Cheev start ***\n");
            exit(1);
        }
    }
}

void Cheev_test(void) {
    Cheev_test2("V", "U");
    Cheev_test2("V", "L");
    Cheev_test2("N", "U");
    Cheev_test2("N", "L");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Cheev start ***\n");
    Cheev_test();
    printf("*** Testing Cheev successful ***\n");
    return (0);
}
