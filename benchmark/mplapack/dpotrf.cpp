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

// https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_potrf(int n_i) {
    double adds, muls, flops;
    double n;
    n = (double)n_i;
    muls = (1. / 6.) * n * n * n + 0.5 * n * n + (1. / 3.) * n;
    adds = (1. / 6.) * n * n * n - (1. / 6.) * n;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    char uplo = 'u';
    int STEP = 1, TOTALSTEPS = 400, n = 1;

    double mtemp, dummy;
    double elapsedtime;
    int lda, info;
    int i, j, k, p;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                n = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-U", argv[i]) == 0) {
                uplo = 'u';
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-L", argv[i]) == 0) {
                uplo = 'l';
            }
        }
    }
    for (p = 0; p < TOTALSTEPS; p++) {
        lda = n;
        double *a = new double[lda * n];
        double *a_ref = new double[lda * n];
        for (i = 0; i < lda * n; i++) {
            a[i] = randomnumber(dummy);
        }
        // Positive semidefinite matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mtemp = 0.0;
                for (k = 0; k < n; k++) {
                    mtemp = mtemp + a[i + k * lda] * a[j + k * lda];
                }
                a_ref[i + j * lda] = mtemp;
            }
        }
        for (i = 0; i < lda * n; i++) {
            a[i] = a_ref[i];
        }
        auto t1 = Clock::now();
        dpotrf_f77(&uplo, &n, a, &lda, &info);
        auto t2 = Clock::now();
        elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        printf("    n     MFLOPS   uplo\n");
        printf("%5d %10.3f      %c\n", (int)n, flops_potrf(n) / elapsedtime * MFLOPS, uplo);
        delete[] a_ref;
        delete[] a;
        n = n + STEP;
    }
}
