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

#include <stdio.h>
#include <string.h>
#include <chrono>
#include <dlfcn.h>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_benchmark.h>

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_gemm(mplapackint k_i, mplapackint m_i, mplapackint n_i) {
    double adds, muls, flops;
    double k, m, n;
    m = (double)m_i;
    n = (double)n_i;
    k = (double)k_i;
    muls = m * (k + 2) * n;
    adds = m * k * n;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    REAL alpha, beta, dummy;
    REAL *dummywork;
    double elapsedtime;
    char transa, transb, normtype;
    int N0, M0, K0, STEPN = 3, STEPM = 3, STEPK = 3, LOOPS = 3, TOTALSTEPS = 333;
    int lda, ldb, ldc;
    int i, m, n, k, ka, kb, p;
    int check_flag = 1;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    ___MPLAPACK_INITIALIZE___

    const char mpblas_sym[] = SYMBOL_GCC_RGEMM;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mpblas_ref)(const char *, const char *, mplapackint, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint);
    void (*raxpy_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    // initialization
    N0 = M0 = K0 = 1;
    transa = transb = 'n';
    normtype = 'm';
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-M", argv[i]) == 0) {
                M0 = atoi(argv[++i]);
            } else if (strcmp("-K", argv[i]) == 0) {
                K0 = atoi(argv[++i]);
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
            }
        }
    }

    if (check_flag) {
        handle = dlopen(MPBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) {
            printf("dlopen: %s\n", dlerror());
            return 1;
        }
        mpblas_ref = (void (*)(const char *, const char *, mplapackint, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint))dlsym(handle, mpblas_sym);
        if ((error = dlerror()) != NULL) {
            fprintf(stderr, "%s\n", error);
            return 1;
        }
        raxpy_ref = (void (*)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint))dlsym(handle, raxpy_sym);
        if ((error = dlerror()) != NULL) {
            fprintf(stderr, "%s\n", error);
            return 1;
        }
    }

    m = M0;
    n = N0;
    k = K0;
    for (p = 0; p < TOTALSTEPS; p++) {
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

        REAL *a = new REAL[lda * ka];
        REAL *b = new REAL[ldb * kb];
        REAL *c = new REAL[ldc * n];
        REAL *c_ref = new REAL[ldc * n];
        REAL mOne = -1;
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (i = 0; i < lda * ka; i++) {
            a[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldb * kb; i++) {
            b[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldc * n; i++) {
            c[i] = c_ref[i] = randomnumber(dummy);
        }
        if (check_flag) {
            auto t1 = Clock::now();
            Rgemm(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            (*mpblas_ref)(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c_ref, ldc);
            (*raxpy_ref)((mplapackint)(ldc * n), mOne, c, (mplapackint)1, c_ref, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)ldc, (mplapackint)n, c_ref, ldc, dummywork);
            diffr = cast2double(diff);
            printf("    m     n     k       MFLOPS     error   transa   transb\n");
            printf("%5d %5d %5d  %10.3f    %5.2e       %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, diffr, transa, transb);
        } else {
            elapsedtime = 0.0;
            for (int j = 0; j < LOOPS; j++) {
                auto t1 = Clock::now();
                Rgemm(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
                auto t2 = Clock::now();
                elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            }
            elapsedtime = elapsedtime / (double)LOOPS;
            printf("    m     n     k     MFLOPS    transa   transb\n");
            printf("%5d %5d %5d %10.3f         %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, transa, transb);
        }
        delete[] c_ref;
        delete[] c;
        delete[] b;
        delete[] a;
        m = m + STEPM;
        n = n + STEPN;
        k = k + STEPK;
    }
    if (check_flag)
        dlclose(handle);
}
