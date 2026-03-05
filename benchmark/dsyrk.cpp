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
static double flops_syrk(int k_i, int n_i) {
    double n = (double)n_i, k = (double)k_i;
    double muls = k * n * (n + 1) * 0.5 + n * n + n;
    double adds = k * n * (n + 1) * 0.5;
    return muls + adds;
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
    double alpha, beta, dummy;
    char uplo = 'u', trans = 'n';
    int n = 1, k = 1, STEPN = 3, STEPK = 3, LOOPS = 3, TOTALSTEPS = 340;
    int WARMUP = 1;
    int lda, ldc;
    int ka;
    int csv_flag = 0, stats_flag = 0;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    for (int i = 1; i < argc; i++) {
        if      (strcmp("-N",          argv[i]) == 0) n          = atoi(argv[++i]);
        else if (strcmp("-K",          argv[i]) == 0) k          = atoi(argv[++i]);
        else if (strcmp("-STEPN",      argv[i]) == 0) STEPN      = atoi(argv[++i]);
        else if (strcmp("-STEPK",      argv[i]) == 0) STEPK      = atoi(argv[++i]);
        else if (strcmp("-UN",         argv[i]) == 0) { uplo = 'u'; trans = 'n'; }
        else if (strcmp("-UT",         argv[i]) == 0) { uplo = 'u'; trans = 't'; }
        else if (strcmp("-UC",         argv[i]) == 0) { uplo = 'u'; trans = 'c'; }
        else if (strcmp("-LN",         argv[i]) == 0) { uplo = 'l'; trans = 'n'; }
        else if (strcmp("-LT",         argv[i]) == 0) { uplo = 'l'; trans = 't'; }
        else if (strcmp("-LC",         argv[i]) == 0) { uplo = 'l'; trans = 'c'; }
        else if (strcmp("-LOOPS",      argv[i]) == 0) LOOPS      = atoi(argv[++i]);
        else if (strcmp("-WARMUP",     argv[i]) == 0) WARMUP     = atoi(argv[++i]);
        else if (strcmp("-TOTALSTEPS", argv[i]) == 0) TOTALSTEPS = atoi(argv[++i]);
        else if (strcmp("-CSV",        argv[i]) == 0) csv_flag   = 1;
        else if (strcmp("-STATS",      argv[i]) == 0) stats_flag = 1;
    }

    int max_n = n + STEPN * (TOTALSTEPS - 1);
    int max_k = k + STEPK * (TOTALSTEPS - 1);
    size_t max_lda = (size_t)std::max(max_n, max_k);
    double *a = new double[max_lda * (size_t)std::max(max_n, max_k)];
    double *c = new double[(size_t)max_n * max_n];

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) printf("n,k,uplo,trans,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        else            printf("n,k,uplo,trans,loops,mflops\n");
    } else {
        if (stats_flag) printf("    n     k  MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)  uplo trans loops\n");
        else            printf("    n     k      MFLOPS     uplo   trans   loops\n");
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        if (lsame_f77(&trans, "n")) { ka = k; lda = n; } else { ka = n; lda = k; }
        ldc = n;

        alpha = randomnumber(dummy);
        beta  = randomnumber(dummy);
        for (int i = 0; i < lda * ka; i++) a[i] = randomnumber(dummy);
        for (int i = 0; i < ldc * n;  i++) c[i] = randomnumber(dummy);

        for (int w = 0; w < WARMUP; w++)
            dsyrk_f77(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);

        for (int j = 0; j < LOOPS; j++) {
            for (int i = 0; i < ldc * n; i++) c[i] = randomnumber(dummy);
            auto t1 = Clock::now();
            dsyrk_f77(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        double flop = flops_syrk(k, n);
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag)
                printf("%d,%d,%c,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n", n, k, uplo, trans, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            else
                printf("%5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%    %c    %c   %3d\n", n, k, st.mean, st.median, st.min, st.max, cv, uplo, trans, LOOPS);
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag) printf("%d,%d,%c,%c,%d,%.3f\n", n, k, uplo, trans, LOOPS, mflops_val);
            else          printf("%5d %5d %10.3f      %c      %c     %3d\n", n, k, mflops_val, uplo, trans, LOOPS);
        }

        n += STEPN; k += STEPK;
    }

    delete[] c; delete[] a;
    return 0;
}
