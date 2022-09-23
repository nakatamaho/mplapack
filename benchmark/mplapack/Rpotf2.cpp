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
#include <dlfcn.h>
#include <chrono>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_benchmark.h>

#define TOTALSTEPS 1000

int main(int argc, char *argv[]) {
    REAL alpha, beta, mtemp, dummy;
    REAL *dummywork;
    double elapsedtime;
    char uplo, normtype;
    mplapackint N0, STEP, info;
    mplapackint lda;
    int i, j, m, n, k, ka, kb, p, q;
    int check_flag = 1;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    ___MPLAPACK_INITIALIZE___

    const char mplapack_sym[] = SYMBOL_GCC_RPOTF2;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mplapack_ref)(const char *, mplapackint, REAL *, mplapackint, mplapackint *);
    void (*raxpy_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    // initialization
    N0 = 1;
    STEP = 1;
    uplo = 'u';
    normtype = 'm';
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-U", argv[i]) == 0) {
                uplo = 'u';
            } else if (strcmp("-N0", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-L", argv[i]) == 0) {
                uplo = 'l';
            } else if (strcmp("-NOCHECK", argv[i]) == 0) {
                check_flag = 0;
            }
        }
    }

    if (check_flag) {
        handle = dlopen(MPLAPACK_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
        if (!handle) {
            printf("dlopen: %s\n", dlerror());
            return 1;
        }
        mplapack_ref = (void (*)(const char *, mplapackint, REAL *, mplapackint, mplapackint *))dlsym(handle, mplapack_sym);
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

    n = N0;
    for (p = 0; p < TOTALSTEPS; p++) {
        lda = n;
        REAL *A = new REAL[lda * n];
        REAL *Ad = new REAL[lda * n];
        REAL mOne = -1;
        for (i = 0; i < lda * n; i++) {
            A[i] = randomnumber(dummy);
        }
        // Positive semidefinite matrix
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mtemp = 0.0;
                for (k = 0; k < n; k++) {
                    mtemp = mtemp + A[i + k * lda] * A[j + k * lda];
                }
                Ad[i + j * lda] = mtemp;
            }
        }
        for (i = 0; i < lda * n; i++) {
            A[i] = Ad[i];
        }

        if (check_flag) {
            auto t1 = Clock::now();
            Rpotf2(&uplo, n, A, lda, info);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            (*mplapack_ref)(&uplo, n, Ad, lda, &info);
            (*raxpy_ref)((mplapackint)(lda * n), mOne, A, (mplapackint)1, Ad, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)lda, (mplapackint)n, Ad, lda, dummywork);
            diffr = cast2double(diff);
            printf("    n     MFLOPS     error   uplo\n");
            printf("%5d %10.3f %5.2e      %c\n", (int)n, ((double)n * (double)n * (double)n / 3.0) / elapsedtime * MFLOPS, diffr, uplo);
        } else {
            auto t1 = Clock::now();
            Rpotf2(&uplo, n, A, lda, info);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            printf("    n     MFLOPS   uplo\n");
            printf("%5d %10.3f      %c\n", (int)n, ((double)n * (double)n * (double)n / 3.0) / elapsedtime * MFLOPS, uplo);
        }
        delete[] Ad;
        delete[] A;
        n = n + STEP;
    }
    if (check_flag)
        dlclose(handle);
}
