/*
 * Copyright (c) 2008-2026
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

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <vector>
#include <blas.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
static double flops_gemm(int k_i, int m_i, int n_i) {
    double m = (double)m_i, n = (double)n_i, k = (double)k_i;
    return m * (k + 2) * n + m * k * n;
}

struct BenchStats {
    double mean, stddev, min, max, median;
};

static BenchStats compute_stats(std::vector<double> &v) {
    BenchStats s;
    int n = (int)v.size();
    if (n == 0) { s.mean = s.stddev = s.min = s.max = s.median = 0.0; return s; }
    std::sort(v.begin(), v.end());
    s.min = v[0]; s.max = v[n - 1];
    s.median = (n % 2 == 1) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i];
    s.mean = sum / n;
    double var = 0.0;
    for (int i = 0; i < n; i++) { double d = v[i] - s.mean; var += d * d; }
    s.stddev = (n > 1) ? std::sqrt(var / (n - 1)) : 0.0;
    return s;
}

int main(int argc, char *argv[]) {
    int STEPN = 7, STEPM = 7, STEPK = 7, LOOPS = 7, TOTALSTEPS = 428;
    int WARMUP = 1;
    int m = 1, n = 1, k = 1;
    char transa = 'n', transb = 'n';
    int csv_flag = 0, stats_flag = 0;

    double alpha, beta, dummy;
    int ka, kb, lda, ldb, ldc;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    for (int i = 1; i < argc; i++) {
        if      (strcmp("-N",          argv[i]) == 0) n          = atoi(argv[++i]);
        else if (strcmp("-M",          argv[i]) == 0) m          = atoi(argv[++i]);
        else if (strcmp("-K",          argv[i]) == 0) k          = atoi(argv[++i]);
        else if (strcmp("-STEPN",      argv[i]) == 0) STEPN      = atoi(argv[++i]);
        else if (strcmp("-STEPM",      argv[i]) == 0) STEPM      = atoi(argv[++i]);
        else if (strcmp("-STEPK",      argv[i]) == 0) STEPK      = atoi(argv[++i]);
        else if (strcmp("-NN",         argv[i]) == 0) { transa = transb = 'n'; }
        else if (strcmp("-TT",         argv[i]) == 0) { transa = transb = 't'; }
        else if (strcmp("-NT",         argv[i]) == 0) { transa = 'n'; transb = 't'; }
        else if (strcmp("-TN",         argv[i]) == 0) { transa = 't'; transb = 'n'; }
        else if (strcmp("-LOOPS",      argv[i]) == 0) LOOPS      = atoi(argv[++i]);
        else if (strcmp("-WARMUP",     argv[i]) == 0) WARMUP     = atoi(argv[++i]);
        else if (strcmp("-TOTALSTEPS", argv[i]) == 0) TOTALSTEPS = atoi(argv[++i]);
        else if (strcmp("-CSV",        argv[i]) == 0) csv_flag   = 1;
        else if (strcmp("-STATS",      argv[i]) == 0) stats_flag = 1;
    }

    // Pre-allocate at maximum size
    int max_m = m + STEPM * (TOTALSTEPS - 1);
    int max_n = n + STEPN * (TOTALSTEPS - 1);
    int max_k = k + STEPK * (TOTALSTEPS - 1);
    size_t max_lda = (size_t)std::max(max_m, max_k);
    size_t max_ldb = (size_t)std::max(max_k, max_n);
    size_t max_ldc = (size_t)max_m;
    double *a = new double[max_lda * (size_t)std::max(max_m, max_k)];
    double *b = new double[max_ldb * (size_t)std::max(max_k, max_n)];
    double *c = new double[max_ldc * (size_t)max_n];

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) printf("m,n,k,transa,transb,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        else            printf("m,n,k,transa,transb,loops,mflops\n");
    } else {
        if (stats_flag) printf("    m     n     k   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)     transa transb loops\n");
        else            printf("    m     n     k       MFLOPS   transa   transb   loops\n");
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        if (lsame_f77(&transa, "n")) { ka = k; lda = m; } else { ka = m; lda = k; }
        if (lsame_f77(&transb, "n")) { kb = n; ldb = k; } else { kb = k; ldb = n; }
        ldc = m;

        alpha = randomnumber(dummy);
        beta  = randomnumber(dummy);
        for (int i = 0; i < lda * ka; i++) a[i] = randomnumber(dummy);
        for (int i = 0; i < ldb * kb; i++) b[i] = randomnumber(dummy);
        for (int i = 0; i < ldc * n;  i++) c[i] = randomnumber(dummy);

        for (int w = 0; w < WARMUP; w++)
            dgemm_f77(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);

        for (int j = 0; j < LOOPS; j++) {
            for (int i = 0; i < ldc * n; i++) c[i] = randomnumber(dummy);
            auto t1 = Clock::now();
            dgemm_f77(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        double flop = flops_gemm(k, m, n);
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag)
                printf("%d,%d,%d,%c,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n", m, n, k, transa, transb, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            else
                printf("%5d %5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%      %c      %c   %3d\n", m, n, k, st.mean, st.median, st.min, st.max, cv, transa, transb, LOOPS);
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag)
                printf("%d,%d,%d,%c,%c,%d,%.3f\n", m, n, k, transa, transb, LOOPS, mflops_val);
            else
                printf("%5d %5d %5d  %10.3f        %c        %c     %3d\n", m, n, k, mflops_val, transa, transb, LOOPS);
        }

        m += STEPM; n += STEPN; k += STEPK;
    }

    delete[] c; delete[] b; delete[] a;
    return 0;
}
