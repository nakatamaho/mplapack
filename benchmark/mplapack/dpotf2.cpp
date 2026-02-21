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

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <vector>
#include <blas.h>
#include <lapack.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

// https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
static double flops_potrf(int n_i) {
    double n = (double)n_i;
    double muls = (1.0 / 6.0) * n * n * n + 0.5 * n * n + (1.0 / 3.0) * n;
    double adds = (1.0 / 6.0) * n * n * n - (1.0 / 6.0) * n;
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

// Build a positive semidefinite matrix: A <- A * A^T
static void make_psd(double *a, double *a_ref, int n, int lda) {
    double mtemp;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mtemp = 0.0;
            for (int k = 0; k < n; k++)
                mtemp += a[i + k * lda] * a[j + k * lda];
            a_ref[i + j * lda] = mtemp;
        }
    }
    for (int i = 0; i < lda * n; i++) a[i] = a_ref[i];
}

int main(int argc, char *argv[]) {
    char uplo = 'u';
    int STEP = 1, TOTALSTEPS = 400, n = 1;
    int LOOPS = 3, WARMUP = 1;
    int csv_flag = 0, stats_flag = 0;

    double mtemp, dummy;
    int lda, info;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    for (int i = 1; i < argc; i++) {
        if      (strcmp("-N",          argv[i]) == 0) n          = atoi(argv[++i]);
        else if (strcmp("-STEP",       argv[i]) == 0) STEP       = atoi(argv[++i]);
        else if (strcmp("-U",          argv[i]) == 0) uplo       = 'u';
        else if (strcmp("-L",          argv[i]) == 0) uplo       = 'l';
        else if (strcmp("-LOOPS",      argv[i]) == 0) LOOPS      = atoi(argv[++i]);
        else if (strcmp("-WARMUP",     argv[i]) == 0) WARMUP     = atoi(argv[++i]);
        else if (strcmp("-TOTALSTEPS", argv[i]) == 0) TOTALSTEPS = atoi(argv[++i]);
        else if (strcmp("-CSV",        argv[i]) == 0) csv_flag   = 1;
        else if (strcmp("-STATS",      argv[i]) == 0) stats_flag = 1;
    }

    int max_n = n + STEP * (TOTALSTEPS - 1);
    double *a     = new double[(size_t)max_n * max_n];
    double *a_ref = new double[(size_t)max_n * max_n];
    double *a_psd = new double[(size_t)max_n * max_n]; // saved PSD matrix for re-use across LOOPS

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) printf("n,uplo,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        else            printf("n,uplo,loops,mflops\n");
    } else {
        if (stats_flag) printf("    n  MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)  uplo loops\n");
        else            printf("    n       MFLOPS   uplo loops\n");
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        lda = n;

        // Build PSD matrix once per step
        for (int i = 0; i < lda * n; i++) a[i] = randomnumber(dummy);
        make_psd(a, a_ref, n, lda);
        for (int i = 0; i < lda * n; i++) a_psd[i] = a[i]; // save for LOOPS re-use

        // Warmup
        for (int w = 0; w < WARMUP; w++) {
            for (int i = 0; i < lda * n; i++) a[i] = a_psd[i];
            dpotf2_f77(&uplo, &n, a, &lda, &info);
        }

        // Timed runs
        for (int j = 0; j < LOOPS; j++) {
            for (int i = 0; i < lda * n; i++) a[i] = a_psd[i];
            auto t1 = Clock::now();
            dpotf2_f77(&uplo, &n, a, &lda, &info);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        double flop = flops_potrf(n);
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag)
                printf("%d,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n", n, uplo, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            else
                printf("%5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%    %c  %3d\n", n, st.mean, st.median, st.min, st.max, cv, uplo, LOOPS);
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag) printf("%d,%c,%d,%.3f\n", n, uplo, LOOPS, mflops_val);
            else          printf("%5d  %10.3f      %c  %3d\n", n, mflops_val, uplo, LOOPS);
        }

        n += STEP;
    }

    delete[] a_psd;
    delete[] a_ref;
    delete[] a;
    return 0;
}
