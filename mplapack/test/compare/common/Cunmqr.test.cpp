/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cungqr.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N 10
#define MIN_M 1
#define MAX_M 10
#define MIN_K 1
#define MAX_K 10
#define MAX_LDA 10
#define MAX_LDC 10
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Cunmqr_test2(const char *side, const char *trans) {
    int errorflag = FALSE;
    int iter;
    int m, n, k;
    int lda, ldc, minlda;
    REAL_REF diff;
    INTEGER_REF info_ref, worksize_ref, lwork;
    INTEGER info, worksize;
      
    for (m = MIN_M; m <= MAX_M; m++) {
        for (n = MIN_N; n <= m; n++) {
            for (k = MIN_K; k <= n; k++) {
                if (Mlsame(side, "R")) minlda = max(1, n);
                if (Mlsame(side, "L")) minlda = max(1, m);
                for (lda = minlda; lda <= MAX_LDA; lda++) {
                    for (ldc = max(1, m); ldc <= MAX_LDC; ldc++) {
#if defined VERBOSE_TEST
                        printf("# m %d n %d k %d lda %d\n", m, n, k, lda);
#endif
                        INTEGER_REF lwork_ref = veclen(n, 1) * 1024;
                        COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, k)];
                        COMPLEX_REF *C_ref = new COMPLEX_REF[matlen(ldc, n)];
                        COMPLEX_REF *tau_ref = new COMPLEX_REF[veclen(k, 1)];
                        COMPLEX_REF *work_ref = new COMPLEX_REF[lwork_ref];

                        INTEGER lwork = veclen(n, 1) * 1024;
                        COMPLEX *A = new COMPLEX[matlen(lda, k)];
                        COMPLEX *C = new COMPLEX[matlen(ldc, n)];
                        COMPLEX *tau = new COMPLEX[veclen(k, 1)];
                        COMPLEX *work = new COMPLEX[lwork];

                        for (iter = 0; iter < MAX_ITER; iter++) {
                            set_random_vector(A_ref, A, matlen(lda, k));
                            set_random_vector(C_ref, C, matlen(ldc, n));
                            set_random_vector(tau_ref, tau, veclen(k, 1));
                            set_random_vector(work_ref, work, lwork_ref);
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
                            zunmqr_f77(side, trans, &m, &n, &k, A_ref, &lda, tau_ref, C_ref, &ldc, work_ref, &lwork_ref, &info_ref);
#else
                            Cunmqr(side, trans, m, n, k, A_ref, lda, tau_ref, C_ref, &ldc, work_ref, lwork, info_ref);
#endif
                            Cunmqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);

//                            printf("C="); printmat(ldc,n,C,ldc); printf("\n");
//                            printf("Cref="); printmat(ldc,n,C_ref,ldc); printf("\n");
                            diff = infnorm(C_ref, C, matlen(ldc, n), 1);
                            if (diff > EPSILON11) {
                                printf("error in C: ");
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
                        }
                        delete[] tau_ref;
                        delete[] work_ref;
                        delete[] A_ref;
                        delete[] tau;
                        delete[] work;
                        delete[] A;
                    }
                }
            }
        }
    }
    if (errorflag == TRUE) {
        printf("*** Testing Cunmqr failed ***\n");
        exit(1);
    }
}

void Cunmqr_test() {
    Cunmqr_test2("L", "N");
    Cunmqr_test2("L", "C");
    Cunmqr_test2("R", "N");
    Cunmqr_test2("R", "C");
}

int main(int argc, char *argv[]) {
    printf("*** Testing Cunmqr start ***\n");
    Cunmqr_test();
    printf("*** Testing Cunmqr successful ***\n");
    return (0);
}
