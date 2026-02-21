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

#include <stdio.h>
#include <string.h>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <vector>
#include <dlfcn.h>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_benchmark.h>
#include <mplapack_symbol_resolver.h>

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
    mplapackint n = 1, incx = 1, incy = 1;
    int STEP = 97, LOOPS = 3, TOTALSTEPS = 3092;
    int WARMUP = 1;
    int check_flag = 1;
    int csv_flag = 0, stats_flag = 0;
    bool printlib_flag = false;

    REAL dummy, ans, ans_ref;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    ___MPLAPACK_INITIALIZE___

    typedef REAL (*rdot_func_t)(mplapackint, REAL *, mplapackint, REAL *, mplapackint);

    void *handle = nullptr;
    rdot_func_t mpblas_ref = nullptr;
    double diffr = 0.0;

    for (int i = 1; i < argc; i++) {
        if      (strcmp("-N",          argv[i]) == 0) n          = atoi(argv[++i]);
        else if (strcmp("-STEP",       argv[i]) == 0) STEP       = atoi(argv[++i]);
        else if (strcmp("-NOCHECK",    argv[i]) == 0) check_flag = 0;
        else if (strcmp("-LOOPS",      argv[i]) == 0) LOOPS      = atoi(argv[++i]);
        else if (strcmp("-WARMUP",     argv[i]) == 0) WARMUP     = atoi(argv[++i]);
        else if (strcmp("-TOTALSTEPS", argv[i]) == 0) TOTALSTEPS = atoi(argv[++i]);
        else if (strcmp("-CSV",        argv[i]) == 0) csv_flag   = 1;
        else if (strcmp("-STATS",      argv[i]) == 0) stats_flag = 1;
        else if (strcmp("-PRINTLIB",   argv[i]) == 0) printlib_flag = true;
    }

    if (check_flag) {
        handle = dlopen(MPBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) { fprintf(stderr, "dlopen: %s\n", dlerror()); return 1; }
        mpblas_ref = reinterpret_cast<rdot_func_t>(mplapack_resolver::resolve_symbol(handle, "Rdot", printlib_flag));
        if (!mpblas_ref) {
            fprintf(stderr, "Failed to resolve Rdot in %s%s\n", MPBLAS_REF_LIB, DYLIB_SUFFIX);
            dlclose(handle); return 1;
        }
    }

    mplapackint max_n = n + (mplapackint)STEP * (TOTALSTEPS - 1);
    REAL *x = new REAL[max_n];
    REAL *y = new REAL[max_n];

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) {
            if (check_flag) printf("n,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct,error\n");
            else            printf("n,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        } else {
            if (check_flag) printf("n,loops,mflops,error\n");
            else            printf("n,loops,mflops\n");
        }
    } else {
        if (stats_flag) {
            if (check_flag) printf("         n   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)       error   loops\n");
            else            printf("         n   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)     loops\n");
        } else {
            if (check_flag) printf("         n       MFLOPS      error   loops\n");
            else            printf("         n       MFLOPS    loops\n");
        }
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        for (mplapackint i = 0; i < n; i++) { x[i] = randomnumber(dummy); y[i] = randomnumber(dummy); }

        for (int w = 0; w < WARMUP; w++)
            ans = Rdot(n, x, incx, y, incy);

        for (int j = 0; j < LOOPS; j++) {
            auto t1 = Clock::now();
            ans = Rdot(n, x, incx, y, incy);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        if (check_flag) {
            ans_ref = (*mpblas_ref)(n, x, incx, y, incy);
            diffr   = cast2double(ans - ans_ref);
        }

        double flop = 2.0 * (double)n;
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag) {
                if (check_flag) printf("%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2e\n", (int)n, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv, diffr);
                else            printf("%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n",       (int)n, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            } else {
                if (check_flag) printf("%10d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%  %5.2e   %3d\n", (int)n, st.mean, st.median, st.min, st.max, cv, diffr, LOOPS);
                else            printf("%10d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%   %3d\n",         (int)n, st.mean, st.median, st.min, st.max, cv, LOOPS);
            }
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag) {
                if (check_flag) printf("%d,%d,%.3f,%.2e\n", (int)n, LOOPS, mflops_val, diffr);
                else            printf("%d,%d,%.3f\n",       (int)n, LOOPS, mflops_val);
            } else {
                if (check_flag) printf("%10d   %10.3f   %5.2e   %3d\n", (int)n, mflops_val, diffr, LOOPS);
                else            printf("%10d   %10.3f   %3d\n",          (int)n, mflops_val, LOOPS);
            }
        }

        n += STEP;
    }

    delete[] y;
    delete[] x;
    if (check_flag) dlclose(handle);
    return 0;
}
