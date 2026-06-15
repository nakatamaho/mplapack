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
#include <dlfcn.h>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_benchmark.h>
#include <mplapack_symbol_resolver.h>

void Rgemm_binary128_opencl_kernel(bool btransa, bool btransb,
                                   int m, int n, int k,
                                   mplapack_binary128_t alpha,
                                   mplapack_binary128_t *A, int lda,
                                   mplapack_binary128_t *B, int ldb,
                                   mplapack_binary128_t beta,
                                   mplapack_binary128_t *C, int ldc);

static double flops_gemm(mplapackint k_i, mplapackint m_i, mplapackint n_i) {
    double m = (double)m_i;
    double n = (double)n_i;
    double k = (double)k_i;
    double muls = m * (k + 2) * n;
    double adds = m * k * n;
    return muls + adds;
}

int main(int argc, char *argv[]) {
    REAL alpha, beta, dummy;
    REAL dummywork[1];
    double elapsedtime;
    char transa = 'n', transb = 'n', normtype = 'm';
    int m = 1, n = 1, k = 1, STEPN = 7, STEPM = 7, STEPK = 7, LOOPS = 3, TOTALSTEPS = 720;
    int lda, ldb, ldc;
    int ka, kb;
    int check_flag = 1;
    bool printlib_flag = false;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    typedef void (*rgemm_func_t)(const char *, const char *, mplapackint, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint);
    typedef void (*raxpy_func_t)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);

    void *handle = nullptr;
    rgemm_func_t mpblas_ref = nullptr;
    raxpy_func_t raxpy_ref = nullptr;
    REAL diff;
    double diffr = 0.0;

    ___MPLAPACK_INITIALIZE___

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
        } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
            TOTALSTEPS = atoi(argv[++i]);
        } else if (strcmp("-PRINTLIB", argv[i]) == 0) {
            printlib_flag = true;
        }
    }

    if (check_flag) {
        handle = dlopen(MPBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) {
            fprintf(stderr, "dlopen: %s\n", dlerror());
            return 1;
        }
        mpblas_ref = reinterpret_cast<rgemm_func_t>(mplapack_resolver::resolve_symbol(handle, "Rgemm", printlib_flag));
        if (!mpblas_ref) {
            fprintf(stderr, "Failed to resolve Rgemm\n");
            dlclose(handle);
            return 1;
        }
        raxpy_ref = reinterpret_cast<raxpy_func_t>(mplapack_resolver::resolve_symbol(handle, "Raxpy", printlib_flag));
        if (!raxpy_ref) {
            fprintf(stderr, "Failed to resolve Raxpy\n");
            dlclose(handle);
            return 1;
        }
    }

    printf("    m     n     k     MFLOPS");
    if (check_flag)
        printf("      error");
    printf("   transa   transb\n");

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

        REAL *A = new REAL[lda * ka];
        REAL *B = new REAL[ldb * kb];
        REAL *C = new REAL[ldc * n];
        REAL *Cref = check_flag ? new REAL[ldc * n] : nullptr;
        REAL mOne = -1;

        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (int i = 0; i < lda * ka; i++)
            A[i] = randomnumber(dummy);
        for (int i = 0; i < ldb * kb; i++)
            B[i] = randomnumber(dummy);
        for (int i = 0; i < ldc * n; i++) {
            C[i] = randomnumber(dummy);
            if (check_flag)
                Cref[i] = C[i];
        }

        bool btransa = Mlsame(&transa, "t");
        bool btransb = Mlsame(&transb, "t");

        if (check_flag) {
            auto t1 = Clock::now();
            Rgemm_binary128_opencl_kernel(btransa, btransb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
        } else {
            elapsedtime = 0.0;
            for (int j = 0; j < LOOPS; j++) {
                auto t1 = Clock::now();
                Rgemm_binary128_opencl_kernel(btransa, btransb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
                auto t2 = Clock::now();
                elapsedtime += (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            }
            elapsedtime = elapsedtime / (double)LOOPS;
        }

        if (check_flag) {
            (*mpblas_ref)(&transa, &transb, m, n, k, alpha, A, lda, B, ldb, beta, Cref, ldc);
            (*raxpy_ref)((mplapackint)(ldc * n), mOne, C, (mplapackint)1, Cref, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)ldc, (mplapackint)n, Cref, ldc, dummywork);
            diffr = cast2double(diff);
            printf("%5d %5d %5d  %10.3f    %5.2e     %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, diffr, transa, transb);
        } else {
            printf("%5d %5d %5d  %10.3f        %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, transa, transb);
        }

        delete[] Cref;
        delete[] C;
        delete[] B;
        delete[] A;

        m += STEPM;
        n += STEPN;
        k += STEPK;
        fflush(stdout);
    }

    if (check_flag)
        dlclose(handle);

    return 0;
}
