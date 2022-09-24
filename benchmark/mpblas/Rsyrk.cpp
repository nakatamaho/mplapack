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
double flops_syrk(mplapackint k_i, mplapackint n_i) {
    double adds, muls, flops;
    double n, k;
    n = (double)n_i;
    k = (double)k_i;
    muls = k * n * (n + 1) * 0.5 + n * n + n;
    adds = k * n * (n + 1) * 0.5;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    REAL alpha, beta, dummy;
    REAL dummywork[1];
    double elapsedtime;
    char uplo = 'u', trans = 'n', normtype = 'm';
    int n = 1, k = 1, STEPN = 3, STEPK = 3, LOOPS = 3, TOTALSTEPS = 340;
    int lda, ldc;
    int i, ka, p;
    int check_flag = 1;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    ___MPLAPACK_INITIALIZE___

    const char mpblas_sym[] = SYMBOL_GCC_RSYRK;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mpblas_ref)(const char *, const char *, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL, REAL *, mplapackint);
    void (*raxpy_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    // initialization
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                n = atoi(argv[++i]);
            } else if (strcmp("-K", argv[i]) == 0) {
                k = atoi(argv[++i]);
            } else if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPK", argv[i]) == 0) {
                STEPK = atoi(argv[++i]);
            } else if (strcmp("-UN", argv[i]) == 0) {
                uplo = 'u', trans = 'n';
            } else if (strcmp("-UT", argv[i]) == 0) {
                uplo = 'u', trans = 't';
            } else if (strcmp("-UC", argv[i]) == 0) {
                uplo = 'u', trans = 'c';
            } else if (strcmp("-LN", argv[i]) == 0) {
                uplo = 'l', trans = 'n';
            } else if (strcmp("-LT", argv[i]) == 0) {
                uplo = 'l', trans = 't';
            } else if (strcmp("-LC", argv[i]) == 0) {
                uplo = 'l', trans = 'c';
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
        mpblas_ref = (void (*)(const char *, const char *, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL, REAL *, mplapackint))dlsym(handle, mpblas_sym);
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
    for (p = 0; p < TOTALSTEPS; p++) {
        if (Mlsame(&trans, "n")) {
            ka = k;
            lda = n;
        } else {
            ka = n;
            lda = k;
        }
        ldc = n;

        REAL *a = new REAL[lda * ka];
        REAL *c = new REAL[ldc * n];
        REAL *c_ref = new REAL[ldc * n];
        REAL mOne = -1;
        //    alpha = randomnumber (dummy);
        //    beta = randomnumber (dummy);

        alpha = 1.0;
        beta = 0.0;
        for (i = 0; i < lda * ka; i++) {
            a[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldc * n; i++) {
            c[i] = c_ref[i] = randomnumber(dummy);
        }

        if (check_flag) {
            auto t1 = Clock::now();
            Rsyrk(&uplo, &trans, n, k, alpha, a, lda, beta, c, ldc);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            (*mpblas_ref)(&uplo, &trans, n, k, alpha, a, lda, beta, c_ref, ldc);
            (*raxpy_ref)((mplapackint)(ldc * n), mOne, c, (mplapackint)1, c_ref, (mplapackint)1);

            diff = Rlange(&normtype, (mplapackint)ldc, (mplapackint)n, c_ref, ldc, dummywork);
            diffr = cast2double(diff);
            printf("    n     k      MFLOPS       error    uplo    trans\n");
            printf("%5d %5d  %10.3f    %5.2e       %c        %c\n", (int)n, (int)k, flops_syrk(k, n) / elapsedtime * MFLOPS, diffr, uplo, trans);
        } else {
            elapsedtime = 0.0;
            for (int j = 0; j < LOOPS; j++) {
                auto t1 = Clock::now();
                Rsyrk(&uplo, &trans, n, k, alpha, a, lda, beta, c, ldc);
                auto t2 = Clock::now();
                elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            }
            elapsedtime = elapsedtime / (double)LOOPS;
            printf("    n     k      MFLOPS     uplo   trans\n");
            printf("%5d %5d %10.3f      %c      %c\n", (int)n, (int)k, flops_syrk(k, n) / elapsedtime * MFLOPS, uplo, trans);
        }
        delete[] c_ref;
        delete[] c;
        delete[] a;
        n = n + STEPN;
        k = k + STEPK;
    }
    if (check_flag)
        dlclose(handle);
}
