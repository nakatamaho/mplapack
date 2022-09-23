/*
 * Copyright (c) 2008-2022
 *	Nakata, Maho
 * 	All rights reserved.
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
#include <chrono>
#include <blas.h>
#include <lapack.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

#define TOTALSTEPS 1000

int main(int argc, char *argv[]) {
    double alpha, beta, mtemp, dummy;
    double *dummywork;
    double elapsedtime;
    char uplo, normtype;
    int N0, STEP, info;
    int lda;
    int i, j, m, n, k, ka, kb, p, q;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    N0 = 1;
    STEP = 1;
    uplo = 'u';
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-U", argv[i]) == 0) {
                uplo = 'u';
            } else if (strcmp("-N0", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-L", argv[i]) == 0) {
                uplo = 'l';
            }
        }
    }
    n = N0;
    for (p = 0; p < TOTALSTEPS; p++) {
        lda = n;
        double *A = new double[lda * n];
        double *Ad = new double[lda * n];
        double mOne = -1;
        for (i = 0; i < lda * n; i++) {
            A[i] = randomnumber(dummy);
        }
        // Positive semidefinite matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mtemp = 0.0;
                for (k = 0; k < n; k++) {
                    mtemp = mtemp + A[i + k * lda] * A[j + k * lda];
                }
                Ad[i + j * lda] = mtemp;
            }
        }
        for (i = 0; i < lda * n; i++) {
            A[i] = Ad[i];
        }
        auto t1 = Clock::now();
        dpotrf_f77(&uplo, &n, A, &lda, &info);
        auto t2 = Clock::now();
        elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        printf("    n     MFLOPS   uplo\n");
        printf("%5d %10.3f      %c\n", (int)n, ((double)n * (double)n * (double)n / 3.0) / elapsedtime * MFLOPS, uplo);
        delete[] Ad;
        delete[] A;
        n = n + STEP;
    }
}
