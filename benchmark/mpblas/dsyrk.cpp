/*
 * Copyright (c) 2008-2010
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
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <blas.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_syrk(int k_i, int n_i) {
    double adds, muls, flops;
    double n, k;
    n = (double)n_i;
    k = (double)k_i;
    muls = k * n * (n + 1) * 0.5 + n * n + n;
    adds = k * n * (n + 1) * 0.5;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    double alpha, beta, dummy;
    double elapsedtime;
    char uplo = 'u', trans = 'n';
    int n = 1, k = 1, STEPN = 3, STEPK = 3, LOOPS = 3, TOTALSTEPS = 340;
    int lda, ldc;
    int i, ka, p;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                n = atoi(argv[++i]);
            } else if (strcmp("-K", argv[i]) == 0) {
                k = atoi(argv[++i]);
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
            }
        }
    }
    for (p = 0; p < TOTALSTEPS; p++) {
        if (lsame_f77(&trans, "n")) {
            ka = k;
            lda = n;
        } else {
            ka = n;
            lda = k;
        }
        ldc = n;

        double *a = new double[lda * ka];
        double *c = new double[ldc * n];
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (i = 0; i < lda * ka; i++) {
            a[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldc * n; i++) {
            c[i] = randomnumber(dummy);
        }

        elapsedtime = 0.0;
        for (int j = 0; j < LOOPS; j++) {
            auto t1 = Clock::now();
            dsyrk_f77(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
            auto t2 = Clock::now();
            elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }
        elapsedtime = elapsedtime / (double)LOOPS;
        printf("    n     k      MFLOPS     uplo   trans\n");
        printf("%5d %5d %10.3f      %c      %c\n", (int)n, (int)k, flops_syrk(k, n) / elapsedtime * MFLOPS, uplo, trans);
        delete[] c;
        delete[] a;
        n = n + STEPN;
        k = k + STEPK;
    }
}
