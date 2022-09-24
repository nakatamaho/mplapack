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

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_gemv(int m_i, int n_i) {
    double adds, muls, flops;
    double n, m;
    n = (double)n_i;
    m = (double)m_i;
    muls = m * n + 2. * m;
    adds = m * n;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    int k, l, m, n;
    int STEPN = 7, STEPM = 7, N0, M0, LOOPS = 3, TOTALSTEPS = 283;
    int incx = 1, incy = 1;
    double alpha, beta, dummy, *dummywork;
    double mOne = -1;
    double elapsedtime, t1, t2;
    int i, p;
    char trans, normtype;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    N0 = M0 = 1;
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-M", argv[i]) == 0) {
                M0 = atoi(argv[++i]);
            } else if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPM", argv[i]) == 0) {
                STEPM = atoi(argv[++i]);
            } else if (strcmp("-T", argv[i]) == 0) {
                trans = 't';
            } else if (strcmp("-LOOPS", argv[i]) == 0) {
                LOOPS = atoi(argv[++i]);
            } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
                TOTALSTEPS = atoi(argv[++i]);
            }
        }
    }
    n = N0;
    m = M0;
    for (p = 0; p < TOTALSTEPS; p++) {
        if (lsame_f77(&trans, "n")) {
            k = n;
            l = m;
        } else {
            k = m;
            l = n;
        }
        double *x = new double[k];
        double *y = new double[l];
        double *a = new double[n * m];
        for (i = 0; i < k; i++) {
            x[i] = randomnumber(dummy);
        }
        for (i = 0; i < l; i++) {
            y[i] = randomnumber(dummy);
        }
        for (i = 0; i < k * l; i++) {
            a[i] = randomnumber(dummy);
        }
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        elapsedtime = 0.0;
        for (int j = 0; j < LOOPS; j++) {
            auto t1 = Clock::now();
            dgemv_f77(&trans, &m, &n, &alpha, a, &m, x, &incx, &beta, y, &incy);
            auto t2 = Clock::now();
            elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }
        elapsedtime = elapsedtime / (double)LOOPS;
        printf("     m       n      MFLOPS  trans\n");
        printf("%6d  %6d  %10.3f      %c\n", (int)n, (int)m, flops_gemv(m, n) / elapsedtime * MFLOPS, trans);
        delete[] y;
        delete[] x;
        delete[] a;
        n = n + STEPN;
        m = m + STEPM;
    }
}
