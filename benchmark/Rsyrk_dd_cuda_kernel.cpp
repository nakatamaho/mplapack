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

#include <cuda.h>
#include <cuda_runtime.h>

void Rsyrk_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, dd_real alpha, dd_real *Adev, mplapackint lda, dd_real beta, dd_real *Cdev, mplapackint ldc);

static void SetDevice() {
    int gpudevice = 0;
    int device_count = 0;
    int device;

    cudaGetDeviceCount(&device_count);
    printf("device_count : %d\n", device_count);

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, gpudevice);
    printf("device name -> %s\n", deviceProp.name);

    if (deviceProp.warpSize <= 1) {
        printf("warning, CUDA Device Emulation (CPU) detected, exiting\n");
        exit(1);
    }

    cudaError_t cudareturn = cudaSetDevice(gpudevice);
    printf("cudareturn -> %d\n", cudareturn);
    if (cudareturn == cudaErrorInvalidDevice) {
        perror("cudaSetDevice returned cudaErrorInvalidDevice");
    } else {
        cudaGetDevice(&device);
        printf("cudaGetDevice()=%d\n", device);
    }
}

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
static double flops_syrk(mplapackint k_i, mplapackint n_i) {
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
    REAL alpha, beta, dummy;
    REAL dummywork[1];
    double *dummyd;
    char uplo = 'u', trans = 'n', normtype = 'm';
    int STEPN = 1, STEPK = 1, LOOPS = 3, TOTALSTEPS = 100;
    int WARMUP = 1;
    int lda, ldc;
    int n = 1, k = 1, ka;
    int check_flag = 1;
    int csv_flag = 0, stats_flag = 0;
    bool printlib_flag = false;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    typedef void (*rsyrk_func_t)(const char *, const char *, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL, REAL *, mplapackint);
    typedef void (*raxpy_func_t)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);

    void *handle = nullptr;
    rsyrk_func_t mpblas_ref = nullptr;
    raxpy_func_t raxpy_ref  = nullptr;
    REAL diff;
    double diffr = 0.0;

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
        mpblas_ref = reinterpret_cast<rsyrk_func_t>(mplapack_resolver::resolve_symbol(handle, "Rsyrk", printlib_flag));
        if (!mpblas_ref) { fprintf(stderr, "Failed to resolve Rsyrk\n");  dlclose(handle); return 1; }
        raxpy_ref = reinterpret_cast<raxpy_func_t>(mplapack_resolver::resolve_symbol(handle, "Raxpy", printlib_flag));
        if (!raxpy_ref)  { fprintf(stderr, "Failed to resolve Raxpy\n");  dlclose(handle); return 1; }
    }

    SetDevice();
    // dummy allocation for CUDA context initialization
    cudaMalloc((void **)&dummyd, 16);
    cudaFree(dummyd);

    std::vector<double> times(LOOPS);

    if (csv_flag) {
        if (stats_flag) {
            if (check_flag) printf("n,k,uplo,trans,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct,error\n");
            else            printf("n,k,uplo,trans,loops,mflops_mean,mflops_median,mflops_min,mflops_max,mflops_stddev,cv_pct\n");
        } else {
            if (check_flag) printf("n,k,uplo,trans,loops,mflops,error\n");
            else            printf("n,k,uplo,trans,loops,mflops\n");
        }
    } else {
        if (stats_flag) {
            if (check_flag) printf("    n     k  MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)       error   uplo trans loops\n");
            else            printf("    n     k  MFLOPS(mean) MFLOPS(med)  MFLOPS(min)  MFLOPS(max)  cv(%%)     uplo trans loops\n");
        } else {
            if (check_flag) printf("    n     k      MFLOPS       error    uplo    trans   loops\n");
            else            printf("    n     k      MFLOPS       uplo    trans   loops\n");
        }
    }

    for (int p = 0; p < TOTALSTEPS; p++) {
        if (Mlsame(&trans, "n")) { ka = k; lda = n; } else { ka = n; lda = k; }
        ldc = n;

        REAL *A  = new REAL[lda * ka];
        REAL *C  = new REAL[ldc * n];
        REAL *Cd = check_flag ? new REAL[ldc * n] : nullptr;
        REAL mOne = -1;

        alpha = randomnumber(dummy);
        beta  = randomnumber(dummy);
        for (int i = 0; i < lda * ka; i++) A[i] = randomnumber(dummy);
        for (int i = 0; i < ldc * n;  i++) {
            C[i] = randomnumber(dummy);
            if (check_flag) Cd[i] = C[i];
        }

        int size_A = Mlsame(&trans, "n") ? (lda * k - (lda - n)) : (lda * n - (lda - k));
        int size_C = ldc * n - (ldc - n);

        dd_real *Adev, *Cdev;

        // Warmup
        {
            cudaMalloc((void **)&Adev, size_A * sizeof(REAL));
            cudaMalloc((void **)&Cdev, size_C * sizeof(REAL));
            cudaMemcpy(Adev, A, size_A * sizeof(REAL), cudaMemcpyHostToDevice);
            cudaMemcpy(Cdev, C, size_C * sizeof(REAL), cudaMemcpyHostToDevice);
            for (int w = 0; w < WARMUP; w++)
                Rsyrk_cuda(&uplo, &trans, n, k, alpha, Adev, lda, beta, Cdev, ldc);
            cudaFree(Adev);
            cudaFree(Cdev);
        }

        // Timed runs
        for (int j = 0; j < LOOPS; j++) {
            cudaMalloc((void **)&Adev, size_A * sizeof(REAL));
            cudaMalloc((void **)&Cdev, size_C * sizeof(REAL));
            cudaMemcpy(Adev, A, size_A * sizeof(REAL), cudaMemcpyHostToDevice);
            cudaMemcpy(Cdev, C, size_C * sizeof(REAL), cudaMemcpyHostToDevice);

            auto t1 = Clock::now();
            Rsyrk_cuda(&uplo, &trans, n, k, alpha, Adev, lda, beta, Cdev, ldc);
            auto t2 = Clock::now();
            times[j] = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;

            // Copy result back only on last iteration (for check)
            if (j == LOOPS - 1)
                cudaMemcpy(C, Cdev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
            cudaFree(Adev);
            cudaFree(Cdev);
        }

        if (check_flag) {
            (*mpblas_ref)(&uplo, &trans, n, k, alpha, A, lda, beta, Cd, ldc);
            (*raxpy_ref)((mplapackint)(ldc * n), mOne, C, (mplapackint)1, Cd, (mplapackint)1);
            diff  = Rlange(&normtype, (mplapackint)ldc, (mplapackint)n, Cd, ldc, dummywork);
            diffr = cast2double(diff);
        }

        double flop = flops_syrk(k, n);
        if (stats_flag) {
            std::vector<double> mf(LOOPS);
            for (int j = 0; j < LOOPS; j++) mf[j] = flop / times[j] * MFLOPS;
            BenchStats st = compute_stats(mf);
            double cv = (st.mean > 0.0) ? (st.stddev / st.mean * 100.0) : 0.0;
            if (csv_flag) {
                if (check_flag) printf("%d,%d,%c,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2e\n", n, k, uplo, trans, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv, diffr);
                else            printf("%d,%d,%c,%c,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f\n",       n, k, uplo, trans, LOOPS, st.mean, st.median, st.min, st.max, st.stddev, cv);
            } else {
                if (check_flag) printf("%5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%  %5.2e      %c     %c   %3d\n", n, k, st.mean, st.median, st.min, st.max, cv, diffr, uplo, trans, LOOPS);
                else            printf("%5d %5d  %11.3f  %11.3f  %11.3f  %11.3f  %7.2f%%       %c     %c   %3d\n",        n, k, st.mean, st.median, st.min, st.max, cv, uplo, trans, LOOPS);
            }
        } else {
            BenchStats st = compute_stats(times);
            double mflops_val = flop / st.mean * MFLOPS;
            if (csv_flag) {
                if (check_flag) printf("%d,%d,%c,%c,%d,%.3f,%.2e\n", n, k, uplo, trans, LOOPS, mflops_val, diffr);
                else            printf("%d,%d,%c,%c,%d,%.3f\n",       n, k, uplo, trans, LOOPS, mflops_val);
            } else {
                if (check_flag) printf("%5d %5d  %10.3f    %5.2e       %c        %c     %3d\n", n, k, mflops_val, diffr, uplo, trans, LOOPS);
                else            printf("%5d %5d %10.3f      %c    %c     %3d\n",                  n, k, mflops_val, uplo, trans, LOOPS);
            }
        }

        delete[] Cd;
        delete[] C;
        delete[] A;
        n += STEPN; k += STEPK;
    }

    if (check_flag) dlclose(handle);
    return 0;
}
