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

// https://netlib.org/lapack/lawnspdf/lawn41.pdf
static double flops_getrf(mplapackint m_i, mplapackint n_i) {
    double m = (double)m_i, n = (double)n_i;
    double muls = 0.5 * m * n * n - (1.0 / 6.0) * n * n * n + 0.5 * m * n - 0.5 * n * n + (2.0 / 3.0) * n;
    double adds = 0.5 * m * n * n - (1.0 / 6.0) * n * n * n - 0.5 * m * n + (1.0 / 6.0) * n;
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
    mplapackint m = 1, n = 1;
    mplapackint STEPN = 3, STEPM = 3, TOTALSTEPS = 400;
    int LOOPS = 3, WARMUP = 1;
    int check_flag = 1;
    int csv_flag = 0, stats_flag = 0;
    bool printlib_flag = false;

    char normtype = 'm';
    REAL dummy;
    REAL dummywork[1];
    mplapackint info;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    MPLAPACK_INITIALIZE

    typedef void (*rgetrf_func_t)(mplapackint, mplapackint, REAL *, mplapackint, mplapackint *, mplapackint *);
    typedef void (*raxpy_func_t)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);

    void *handle_mplapack = nullptr;
    void *handle_mpblas   = nullptr;
    rgetrf_func_t mplapack_ref = nullptr;
    raxpy_func_t  raxpy_ref    = nullptr;
    REAL diff;
    double diffr = 0.0;

    for (int i = 1; i < argc; i++) {
        if      (strcmp("-STEPN",      argv[i]) == 0) STEPN      = atoi(argv[++i]);
        else if (strcmp("-STEPM",      argv[i]) == 0) STEPM      = atoi(argv[++i]);
        else if (strcmp("-N",          argv[i]) == 0) n          = atoi(argv[++i]);
        else if (strcmp("-M",          argv[i]) == 0) m          = atoi(argv[++i]);
        else if (strcmp("-NOCHECK",    argv[i]) == 0) check_flag = 0;
        else if (strcmp("-LOOPS",      argv[i]) == 0) LOOPS      = atoi(argv[++i]);
        else if (strcmp("-WARMUP",     argv[i]) == 0) WARMUP     = atoi(argv[++i]);
        else if (strcmp("-TOTALSTEPS", argv[i]) == 0) TOTALSTEPS = atoi(argv[++i]);
        else if (strcmp("-CSV",        argv[i]) == 0) csv_flag   = 1;
        else if (strcmp("-STATS",      argv[i]) == 0) stats_flag = 1;
        else if (strcmp("-PRINTLIB",   argv[i]) == 0) printlib_flag = true;
    }

    if (check_flag) {
        handle_mplapack = dlopen(MPLAPACK_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle_mplapack) { fprintf(stderr, "dlopen: %s\n", dlerror()); return 1; }
        mplapack_ref = reinterpret_cast<rgetrf_func_t>(mplapack_resolver::resolve_symbol(handle_mplapack, "Rgetrf", printlib_flag));
        if (!mplapack_ref) {
            fprintf(stderr, "Failed to resolve Rgetrf in %s%s\n", MPLAPACK_REF_LIB, DYLIB_SUFFIX);
            dlclose(handle_mplapack); return 1;
        }

        handle_mpblas = dlopen(MPBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle_mpblas) { fprintf(stderr, "dlopen: %s\n", dlerror()); dlclose(handle_mplapack); return 1; }
        raxpy_ref = reinterpret_cast<raxpy_func_t>(mplapack_resolver::resolve_symbol(handle_mpblas, "Raxpy", printlib_flag));
        if (!raxpy_ref) {
            fprintf(stderr, "Failed to resolve Raxpy in %s%s\n", MPBLAS_REF_LIB, DYLIB_SUFFIX);
            dlclose(handle_mpblas); dlclose(handle_mplapack); return 1;
        }
    }

    mplapackint max_m   = m + STEPM * (TOTALSTEPS - 1);
    mplapackint max_n   = n + STEPN * (TOTALSTEPS - 1);
    mplapackint max_lda = max_m;
    REAL *a        = new REAL[(size_t)max_lda * max_n];
    REAL *a_ref    = check_flag ? new REAL[(size_t)max_lda * max_n] : nullptr;
    mplapackint *ipiv     = new mplapackint[std::min(max_m, max_n)];
    mplapackint *ipiv_ref = check_flag ? new mplapackint[std::min(max_m, max_n)] : nullptr;
    REAL mOne = -1;

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) {
            if (check_flag) printf("m,n,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct,error\n");
            else            printf("m,n,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        } else {
            if (check_flag) printf("m,n,loops,mflops,error\n");
            else            printf("m,n,loops,mflops\n");
        }
    } else {
        if (stats_flag) {
            if (check_flag) printf("    m     n   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)       error   loops\n");
            else            printf("    m     n   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)     loops\n");
        } else {
            if (check_flag) printf("    m     n       MFLOPS      error   loops\n");
            else            printf("    m     n       MFLOPS    loops\n");
        }
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        mplapackint lda = m;

        for (mplapackint i = 0; i < lda * n; i++) {
            a[i] = randomnumber(dummy);
            if (check_flag) a_ref[i] = a[i];
        }

        // Warmup
        for (int w = 0; w < WARMUP; w++) {
            for (mplapackint i = 0; i < lda * n; i++) a[i] = randomnumber(dummy);
            Rgetrf(m, n, a, lda, ipiv, info);
        }

        // Timed runs
        for (int j = 0; j < LOOPS; j++) {
            for (mplapackint i = 0; i < lda * n; i++) {
                a[i] = randomnumber(dummy);
                if (check_flag) a_ref[i] = a[i];
            }
            auto t1 = Clock::now();
            Rgetrf(m, n, a, lda, ipiv, info);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        // Accuracy check using the last iteration's data
        if (check_flag) {
            (*mplapack_ref)(m, n, a_ref, lda, ipiv_ref, &info);
            (*raxpy_ref)((mplapackint)(lda * n), mOne, a, (mplapackint)1, a_ref, (mplapackint)1);
            diff  = Rlange(&normtype, lda, n, a_ref, lda, dummywork);
            diffr = cast2double(diff);
        }

        double flop = flops_getrf(m, n);
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag) {
                if (check_flag) printf("%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2e\n", (int)m, (int)n, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv, diffr);
                else            printf("%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n",       (int)m, (int)n, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            } else {
                if (check_flag) printf("%5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%  %5.2e   %3d\n", (int)m, (int)n, st.mean, st.median, st.min, st.max, cv, diffr, LOOPS);
                else            printf("%5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%   %3d\n",         (int)m, (int)n, st.mean, st.median, st.min, st.max, cv, LOOPS);
            }
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag) {
                if (check_flag) printf("%d,%d,%d,%.3f,%.2e\n", (int)m, (int)n, LOOPS, mflops_val, diffr);
                else            printf("%d,%d,%d,%.3f\n",       (int)m, (int)n, LOOPS, mflops_val);
            } else {
                if (check_flag) printf("%5d %5d  %10.3f   %5.2e   %3d\n", (int)m, (int)n, mflops_val, diffr, LOOPS);
                else            printf("%5d %5d  %10.3f   %3d\n",          (int)m, (int)n, mflops_val, LOOPS);
            }
        }

        fflush(stdout);
        n += STEPN; m += STEPM;
    }

    delete[] ipiv_ref;
    delete[] ipiv;
    delete[] a_ref;
    delete[] a;
    if (check_flag) { dlclose(handle_mpblas); dlclose(handle_mplapack); }
    return 0;
}
