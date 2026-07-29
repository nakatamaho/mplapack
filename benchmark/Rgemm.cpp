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

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_gemm(mplapackint k_i, mplapackint m_i, mplapackint n_i) {
    double m = (double)m_i;
    double n = (double)n_i;
    double k = (double)k_i;
    double muls = m * (k + 2) * n;
    double adds = m * k * n;
    return muls + adds;
}

// Compute statistics from a vector of samples
struct BenchStats {
    double mean;
    double stddev;
    double min;
    double max;
    double median;
};

static BenchStats compute_stats(std::vector<double> &v) {
    BenchStats s;
    int n = (int)v.size();
    if (n == 0) {
        s.mean = s.stddev = s.min = s.max = s.median = 0.0;
        return s;
    }

    std::sort(v.begin(), v.end());
    s.min = v[0];
    s.max = v[n - 1];
    s.median = (n % 2 == 1) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;

    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += v[i];
    s.mean = sum / n;

    double var = 0.0;
    for (int i = 0; i < n; i++) {
        double d = v[i] - s.mean;
        var += d * d;
    }
    s.stddev = (n > 1) ? std::sqrt(var / (n - 1)) : 0.0;

    return s;
}

int main(int argc, char *argv[]) {
    int STEPN = 3, STEPM = 3, STEPK = 3, LOOPS = 3, TOTALSTEPS = 333;
    int WARMUP = 1;
    int m = 1, n = 1, k = 1;
    char transa = 'n', transb = 'n', normtype = 'm';
    int check_flag = 1;
    int csv_flag = 0;
    int stats_flag = 0;
    bool printlib_flag = false;

    REAL alpha, beta, dummy, dummywork[1];
    int ka, kb, lda, ldb, ldc;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    MPLAPACK_INITIALIZE

    // Typedefs for reference function pointers
    typedef void (*rgemm_func_t)(const char *, const char *, mplapackint, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint);
    typedef void (*raxpy_func_t)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);

    void *handle = nullptr;
    rgemm_func_t mpblas_ref = nullptr;
    raxpy_func_t raxpy_ref = nullptr;
    REAL diff;
    double diffr;

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp("-N", argv[i]) == 0) {
            n = atoi(argv[++i]);
        } else if (strcmp("-M", argv[i]) == 0) {
            m = atoi(argv[++i]);
        } else if (strcmp("-K", argv[i]) == 0) {
            k = atoi(argv[++i]);
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
        } else if (strcmp("-NOCHECK", argv[i]) == 0) {
            check_flag = 0;
        } else if (strcmp("-LOOPS", argv[i]) == 0) {
            LOOPS = atoi(argv[++i]);
        } else if (strcmp("-WARMUP", argv[i]) == 0) {
            WARMUP = atoi(argv[++i]);
        } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
            TOTALSTEPS = atoi(argv[++i]);
        } else if (strcmp("-CSV", argv[i]) == 0) {
            csv_flag = 1;
        } else if (strcmp("-STATS", argv[i]) == 0) {
            stats_flag = 1;
        } else if (strcmp("-PRINTLIB", argv[i]) == 0) {
            printlib_flag = true;
        }
    }

    // Load reference library and resolve symbols dynamically
    if (check_flag) {
        handle = dlopen(MPBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) {
            fprintf(stderr, "dlopen: %s\n", dlerror());
            return 1;
        }
        mpblas_ref = reinterpret_cast<rgemm_func_t>(mplapack_resolver::resolve_symbol(handle, "Rgemm", printlib_flag));
        if (!mpblas_ref) {
            fprintf(stderr, "Failed to resolve Rgemm in %s%s\n", MPBLAS_REF_LIB, DYLIB_SUFFIX);
            dlclose(handle);
            return 1;
        }

        raxpy_ref = reinterpret_cast<raxpy_func_t>(mplapack_resolver::resolve_symbol(handle, "Raxpy", printlib_flag));
        if (!raxpy_ref) {
            fprintf(stderr, "Failed to resolve Raxpy in %s%s\n", MPBLAS_REF_LIB, DYLIB_SUFFIX);
            dlclose(handle);
            return 1;
        }
    }

    // Pre-allocate matrices at maximum size
    int max_m = m + STEPM * (TOTALSTEPS - 1);
    int max_n = n + STEPN * (TOTALSTEPS - 1);
    int max_k = k + STEPK * (TOTALSTEPS - 1);

    size_t max_lda = (size_t)std::max(max_m, max_k);
    size_t max_ka = (size_t)std::max(max_m, max_k);
    size_t max_ldb = (size_t)std::max(max_k, max_n);
    size_t max_kb = (size_t)std::max(max_k, max_n);
    size_t max_ldc = (size_t)max_m;

    size_t size_a = max_lda * max_ka;
    size_t size_b = max_ldb * max_kb;
    size_t size_c = max_ldc * (size_t)max_n;

    REAL *a = new REAL[size_a];
    REAL *b = new REAL[size_b];
    REAL *c = new REAL[size_c];
    REAL *c_ref = check_flag ? new REAL[size_c] : nullptr;
    REAL mOne = -1;

    // Per-loop timing samples
    std::vector<double> times(LOOPS);

    // Print header once
    if (csv_flag) {
        if (stats_flag) {
            if (check_flag)
                printf("m,n,k,transa,transb,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct,error\n");
            else
                printf("m,n,k,transa,transb,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        } else {
            if (check_flag)
                printf("m,n,k,transa,transb,loops,mflops,error\n");
            else
                printf("m,n,k,transa,transb,loops,mflops\n");
        }
    } else {
        if (stats_flag) {
            if (check_flag)
                printf("    m     n     k   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)       error   transa transb loops\n");
            else
                printf("    m     n     k   MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)     transa transb loops\n");
        } else {
            if (check_flag)
                printf("    m     n     k       MFLOPS     error   transa   transb   loops\n");
            else
                printf("    m     n     k       MFLOPS   transa   transb   loops\n");
        }
    }

    // Main benchmark loop
    for (int p = 0; p < TOTALSTEPS; p++) {
        if (Mlsame(&transa, "n")) {
            ka = k;
            lda = m;
        } else {
            ka = m;
            lda = k;
        }
        if (Mlsame(&transb, "n")) {
            kb = n;
            ldb = k;
        } else {
            kb = k;
            ldb = n;
        }
        ldc = m;

        int cn = ldc * n;

        // Fill A and B with random data (once per step)
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (int i = 0; i < lda * ka; i++) {
            a[i] = randomnumber(dummy);
        }
        for (int i = 0; i < ldb * kb; i++) {
            b[i] = randomnumber(dummy);
        }

        // Warmup
        for (int i = 0; i < cn; i++) {
            c[i] = randomnumber(dummy);
        }
        for (int w = 0; w < WARMUP; w++) {
            Rgemm(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        }

        // Timed runs: measure each iteration for per-sample statistics
        for (int j = 0; j < LOOPS; j++) {
            for (int i = 0; i < cn; i++) {
                c[i] = randomnumber(dummy);
                if (check_flag)
                    c_ref[i] = c[i];
            }

            auto t1 = Clock::now();
            Rgemm(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        }

        // Accuracy check using the last iteration
        if (check_flag) {
            (*mpblas_ref)(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c_ref, ldc);
            (*raxpy_ref)((mplapackint)cn, mOne, c, (mplapackint)1, c_ref, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)ldc, (mplapackint)n, c_ref, ldc, dummywork);
            diffr = cast2double(diff);
        }

        double flop = flops_gemm(k, m, n);

        if (stats_flag) {
            // Compute MFLOPS for each sample, then get statistics
            std::vector<double> mflops_v(LOOPS);
            for (int j = 0; j < LOOPS; j++) {
                mflops_v[j] = flop / times[j] * MFLOPS;
            }
            BenchStats st = compute_stats(mflops_v);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;

            if (csv_flag) {
                if (check_flag)
                    printf("%d,%d,%d,%c,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2e\n", m, n, k, transa, transb, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv, diffr);
                else
                    printf("%d,%d,%d,%c,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n", m, n, k, transa, transb, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            } else {
                if (check_flag)
                    printf("%5d %5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%  %5.2e      %c      %c   %3d\n", m, n, k, st.mean, st.median, st.min, st.max, cv, diffr, transa, transb, LOOPS);
                else
                    printf("%5d %5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%       %c      %c   %3d\n", m, n, k, st.mean, st.median, st.min, st.max, cv, transa, transb, LOOPS);
            }
        } else {
            // Simple output: use mean time
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;

            if (csv_flag) {
                if (check_flag)
                    printf("%d,%d,%d,%c,%c,%d,%.3f,%.2e\n", m, n, k, transa, transb, LOOPS, mflops_val, diffr);
                else
                    printf("%d,%d,%d,%c,%c,%d,%.3f\n", m, n, k, transa, transb, LOOPS, mflops_val);
            } else {
                if (check_flag)
                    printf("%5d %5d %5d  %10.3f    %5.2e       %c        %c     %3d\n", m, n, k, mflops_val, diffr, transa, transb, LOOPS);
                else
                    printf("%5d %5d %5d  %10.3f         %c        %c     %3d\n", m, n, k, mflops_val, transa, transb, LOOPS);
            }
        }

        m += STEPM;
        n += STEPN;
        k += STEPK;
    }

    // Cleanup
    delete[] c_ref;
    delete[] c;
    delete[] b;
    delete[] a;

    if (check_flag)
        dlclose(handle);

    return 0;
}
