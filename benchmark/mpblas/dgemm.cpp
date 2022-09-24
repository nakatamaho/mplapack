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
#include <chrono>
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
    double elapsedtime;
    char transa, transb, normtype;
    int N0, M0, K0, STEPN = 7, STEPM = 7, STEPK = 7, LOOP = 7, TOTALSTEPS = 428;
    int lda, ldb, ldc;
    int i, j, m, n, k, ka, kb, p, q;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    N0 = M0 = K0 = 1;
    STEPM = STEPN = STEPK = 1;
    transa = transb = 'n';
    normtype = 'm';
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-M", argv[i]) == 0) {
                M0 = atoi(argv[++i]);
            } else if (strcmp("-K", argv[i]) == 0) {
                K0 = atoi(argv[++i]);
            } else if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPM", argv[i]) == 0) {
                STEPM = atoi(argv[++i]);
            } else if (strcmp("-STEPK", argv[i]) == 0) {
                STEPK = atoi(argv[++i]);
            } else if (strcmp("-NN", argv[i]) == 0) {
                transa = transb = 'n';
            } else if (strcmp("-TT", argv[i]) == 0) {
                transa = transb = 't';
            } else if (strcmp("-NT", argv[i]) == 0) {
                transa = 'n';
                transb = 't';
            } else if (strcmp("-TN", argv[i]) == 0) {
                transa = 't';
                transb = 'n';
            } else if (strcmp("-LOOPS", argv[i]) == 0) {
                LOOPS = atoi(argv[++i]);
            } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
                TOTALSTEPS = atoi(argv[++i]);
            }
        }
    }

    m = M0;
    n = N0;
    k = K0;
    for (p = 0; p < TOTALSTEPS; p++) {
        if (lsame_f77(&transa, "n")) {
            ka = k;
            lda = m;
        } else {
            ka = m;
            lda = k;
        }
        if (lsame_f77(&transb, "n")) {
            kb = n;
            ldb = k;
        } else {
            kb = k;
            ldb = n;
        }
        ldc = m;

        double *a = new double[lda * ka];
        double *b = new double[ldb * kb];
        double *c = new double[ldc * n];
        double mOne = -1;
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (i = 0; i < lda * ka; i++) {
            a[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldb * kb; i++) {
            b[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldc * n; i++) {
            c[i] = randomnumber(dummy);
        }
        auto t1 = Clock::now();
        dgemm_f77(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
        auto t2 = Clock::now();
        elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        printf("    m     n     k     MFLOPS    transa   transb\n");
        printf("%5d %5d %5d %10.3f         %c        %c\n", (int)m, (int)n, (int)k, (2.0 * (double)m * (double)n * (double)k + 2.0 * (double)m * (double)n) / elapsedtime * MFLOPS, transa, transb);
        delete[] c;
        delete[] b;
        delete[] a;
        m = m + STEPM;
        n = n + STEPN;
        k = k + STEPK;
    }
}
