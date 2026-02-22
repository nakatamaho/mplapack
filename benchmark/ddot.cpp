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
    int n = 1, incx = 1, incy = 1, STEP = 97, LOOPS = 3, TOTALSTEPS = 3092;
    int WARMUP = 1;
    int csv_flag = 0, stats_flag = 0;

    double dummy;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    for (int i = 1; i < argc; i++) {
        if      (strcmp("-N",          argv[i]) == 0) n          = atoi(argv[++i]);
        else if (strcmp("-STEP",       argv[i]) == 0) STEP       = atoi(argv[++i]);
        else if (strcmp("-LOOPS",      argv[i]) == 0) LOOPS      = atoi(argv[++i]);
        else if (strcmp("-WARMUP",     argv[i]) == 0) WARMUP     = atoi(argv[++i]);
        else if (strcmp("-TOTALSTEPS", argv[i]) == 0) TOTALSTEPS = atoi(argv[++i]);
        else if (strcmp("-CSV",        argv[i]) == 0) csv_flag   = 1;
        else if (strcmp("-STATS",      argv[i]) == 0) stats_flag = 1;
    }

    int max_n = n + STEP * (TOTALSTEPS - 1);
    double *x = new double[max_n];
    double *y = new double[max_n];

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) printf("n,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        else            printf("n,loops,mflops\n");
    } else {
        if (stats_flag) printf("         n   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)     loops\n");
        else            printf("         n       MFLOPS    loops\n");
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        for (int i = 0; i < n; i++) { x[i] = randomnumber(dummy); y[i] = randomnumber(dummy); }

        for (int w = 0; w < WARMUP; w++)
            ddot_f77(&n, x, &incx, y, &incy);

        for (int j = 0; j < LOOPS; j++) {
            auto t1 = Clock::now();
            ddot_f77(&n, x, &incx, y, &incy);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        double flop = 2.0 * (double)n;
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag)
                printf("%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n", n, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            else
                printf("%10d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%   %3d\n", n, st.mean, st.median, st.min, st.max, cv, LOOPS);
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag) printf("%d,%d,%.3f\n", n, LOOPS, mflops_val);
            else          printf("%10d   %10.3f   %3d\n", n, mflops_val, LOOPS);
        }
        n += STEP;
    }

    delete[] y;
    delete[] x;
    return 0;
}
