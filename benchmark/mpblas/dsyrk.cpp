/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: dgemm.cpp,v 1.5 2010/08/19 01:29:39 nakatamaho Exp $
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

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <blas.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

#define TOTALSTEPS 1000

int main(int argc, char *argv[]) {
    double alpha, beta, dummy;
    double *dummywork;
    double elapsedtime, t1, t2;
    char uplo, trans, normtype;
    int N0, K0, STEPN, STEPK;
    int lda, ldc;
    int i, j, n, k, ka, kb, p, q;
    int check_flag = 1;

    // initialization
    N0 = K0 = 1;
    STEPN = STEPK = 1;
    uplo = 'u';
    trans = 'n';
    normtype = 'm';
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-K", argv[i]) == 0) {
                K0 = atoi(argv[++i]);
            } else if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPK", argv[i]) == 0) {
                STEPK = atoi(argv[++i]);
            } else if (strcmp("-UN", argv[i]) == 0) {
                uplo = 'u', trans = 'n';
            } else if (strcmp("-UT", argv[i]) == 0) {
                uplo = 'u', trans = 't';
            } else if (strcmp("-UC", argv[i]) == 0) {
                uplo = 'u', trans = 'c';
            } else if (strcmp("-LN", argv[i]) == 0) {
                uplo = 'l', trans = 'n';
            } else if (strcmp("-LT", argv[i]) == 0) {
                uplo = 'l', trans = 't';
            } else if (strcmp("-LC", argv[i]) == 0) {
                uplo = 'l', trans = 'c';
            } else if (strcmp("-NOCHECK", argv[i]) == 0) {
                check_flag = 0;
            }
        }
    }

    n = N0;
    k = K0;
    for (p = 0; p < TOTALSTEPS; p++) {
        if (lsame_f77(&trans, "n")) {
            ka = k;
            lda = n;
        } else {
            ka = n;
            lda = k;
        }
        ldc = n;

        double *A = new double[lda * ka];
        double *C = new double[ldc * n];
        double mOne = -1;
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (i = 0; i < lda * ka; i++) {
            A[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldc * n; i++) {
            C[i] = randomnumber(dummy);
        }
        t1 = gettime();
        dsyrk_f77(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
        t2 = gettime();
        elapsedtime = (t2 - t1);
        printf("    n     k      MFLOPS       uplo    trans\n");
        // 2n^2k+2n^2 flops are needed
        printf("%5d %5d  %10.3f     %c        %c\n", (int)n, (int)k, (2.0 * (double)n * (double)n * (double)k + 2.0 * (double)n * (double)n) / elapsedtime * MFLOPS, uplo, trans);
        delete[] C;
        delete[] A;
        n = n + STEPN;
        k = k + STEPK;
    }
}
