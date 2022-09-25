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

// https://netlib.org/lapack/lawnspdf/lawn41.pdf
double flops_getrf(mplapackint m_i, mplapackint n_i) {
    double adds, muls, flops;
    double m, n;
    m = (double)m_i;
    n = (double)n_i;
    muls = 0.5 * m * n * n - (1. / 6.) * n * n * n + 0.5 * m * n - 0.5 * n * n + (2. / 3.) * n;
    adds = 0.5 * m * n * n - (1. / 6.) * n * n * n - 0.5 * m * n + (1. / 6.) * n;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    int m = 1, n = 1;
    mplapackint STEPN = 3, STEPM = 3, TOTALSTEPS = 400;

    char normtype = 'm';
    int check_flag = 1;

    REAL dummy;
    REAL *dummywork = new REAL[1];
    double elapsedtime;
    mplapackint lda, info;
    int i, p;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    ___MPLAPACK_INITIALIZE___

    const char mplapack_sym[] = SYMBOL_GCC_RGETRF;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mplapack_ref)(mplapackint, mplapackint, REAL *, mplapackint, mplapackint *, mplapackint *);
    void (*raxpy_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    // initialization
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPM", argv[i]) == 0) {
                STEPM = atoi(argv[++i]);
            } else if (strcmp("-N", argv[i]) == 0) {
                n = atoi(argv[++i]);
            } else if (strcmp("-M", argv[i]) == 0) {
                m = atoi(argv[++i]);
            } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
                TOTALSTEPS = atoi(argv[++i]);
            }
        }
    }

    if (check_flag) {
        handle = dlopen(MPLAPACK_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) {
            printf("dlopen: %s\n", dlerror());
            return 1;
        }
        mplapack_ref = (void (*)(mplapackint, mplapackint, REAL *, mplapackint, mplapackint *, mplapackint *))dlsym(handle, mplapack_sym);
        if ((error = dlerror()) != NULL) {
            fprintf(stderr, "%s\n", error);
            return 1;
        }

        handle = dlopen(MPBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) {
            printf("dlopen: %s\n", dlerror());
            return 1;
        }
        raxpy_ref = (void (*)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint))dlsym(handle, raxpy_sym);
        if ((error = dlerror()) != NULL) {
            fprintf(stderr, "%s\n", error);
            return 1;
        }
    }
    for (p = 0; p < TOTALSTEPS; p++) {
        lda = m;
        REAL *a = new REAL[lda * n];
        REAL *a_ref = new REAL[lda * n];
        mplapackint *ipiv = new mplapackint[min(m, n)];
        mplapackint *ipiv_ref = new mplapackint[min(m, n)];
        REAL mOne = -1;
        for (i = 0; i < lda * n; i++) {
            a[i] = a_ref[i] = randomnumber(dummy);
        }

        if (check_flag) {
            auto t1 = Clock::now();
            Rgetrf(m, n, a, lda, ipiv, info);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            (*mplapack_ref)(m, n, a_ref, lda, ipiv_ref, &info);
            (*raxpy_ref)((mplapackint)(lda * n), mOne, a, (mplapackint)1, a_ref, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)lda, (mplapackint)n, a_ref, lda, dummywork);
            diffr = cast2double(diff);
            printf("    n     m     MFLOPS   error\n");
            printf("%5d %5d %10.3f %5.2e\n", (int)n, (int)m, flops_getrf(m, n) / elapsedtime * MFLOPS, diffr);
        } else {
            auto t1 = Clock::now();
            Rgetrf(m, n, a, lda, ipiv, info);
            auto t2 = Clock::now();
            elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            printf("    n     m     MFLOPS\n");
            printf("%5d %5d %10.3f\n", (int)n, (int)m, flops_getrf(m, n) / elapsedtime * MFLOPS);
        }
        delete[] ipiv_ref;
        delete[] ipiv;
        delete[] a_ref;
        delete[] a;
        n = n + STEPN;
        m = m + STEPM;
    }
    if (check_flag)
        dlclose(handle);
}
