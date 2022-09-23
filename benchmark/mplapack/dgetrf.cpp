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
#include <string.h>
#include <chrono>
#include <blas.h>
#include <lapack.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

// https://netlib.org/lapack/lawnspdf/lawn41.pdf
double flops_getrf(int m_i, int n_i) {
    double adds, muls, flops;
    double m, n;
    m = (double)m_i;
    n = (double)n_i;
    muls = 0.5 * m * n * n - (1. / 6.) * n * n * n + 0.5 * m * n - 0.5 * n * n + (2. / 3.) * n;
    adds = 0.5 * m * n * n - (1. / 6.) * n * n * n - 0.5 * m * n + (1. / 6.) * n;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    double alpha, beta, dummy;
    double *dummywork;
    double elapsedtime;
    char uplo, normtype;
    int N0, M0, K0, STEPN = 3, STEPM = 3, STEPK = 3, TOTALSTEPS = 400;
    int lda, ldb, ldc, info;
    int i, j, m, n, k, ka, kb, p, q;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    N0 = 1;
    M0 = 1;
    STEPN = 1;
    STEPM = 1;
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPM", argv[i]) == 0) {
                STEPM = atoi(argv[++i]);
            } else if (strcmp("-N0", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-M0", argv[i]) == 0) {
                M0 = atoi(argv[++i]);
            } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
                TOTALSTEPS = atoi(argv[++i]);
            }
        }
    }

    n = N0;
    m = M0;
    for (p = 0; p < TOTALSTEPS; p++) {
        lda = m;
        double *a = new double[lda * n];
        int *ipiv = new int[std::min(m, n)];
        for (i = 0; i < lda * n; i++) {
            a[i] = randomnumber(dummy);
        }
        auto t1 = Clock::now();
        dgetrf_f77(&m, &n, a, &lda, ipiv, &info);
        auto t2 = Clock::now();
        elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        printf("    n     m     MFLOPS\n");
        printf("%5d %5d %10.3f\n", (int)n, (int)m, flops_getrf(m, n) / elapsedtime * MFLOPS);
        delete[] ipiv;
        delete[] a;
        n = n + STEPN;
        m = m + STEPM;
    }
}
