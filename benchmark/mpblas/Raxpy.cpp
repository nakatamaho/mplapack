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

int main(int argc, char *argv[]) {
    mplapackint n = 1;
    mplapackint incx = 1, incy = 1, STEP = 97, LOOPS = 3, TOTALSTEPS = 3092;

    char normtype = 'm';
    int check_flag = 1;

    REAL alpha, dummy, dummywork[1];
    double elapsedtime;
    int i, p;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    ___MPLAPACK_INITIALIZE___

    const char mpblas_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mpblas_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                n = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
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
        mpblas_ref = (void (*)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint))dlsym(handle, mpblas_sym);
        if ((error = dlerror()) != NULL) {
            fprintf(stderr, "%s\n", error);
            return 1;
        }
    }
    for (p = 0; p < TOTALSTEPS; p++) {
        REAL *x = new REAL[n];
        REAL *y = new REAL[n];
        REAL *y_ref = new REAL[n];
        REAL *z = new REAL[n];
        if (check_flag) {
            for (i = 0; i < n; i++) {
                x[i] = randomnumber(dummy);
                y[i] = y_ref[i] = randomnumber(dummy);
            }
            alpha = randomnumber(dummy);
            auto t1 = Clock::now();
            Raxpy(n, alpha, x, incx, y, incy);
            auto t2 = Clock::now();
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            (*mpblas_ref)(n, alpha, x, incx, y_ref, incy);
            for (i = 0; i < n; i++) {
                z[i] = y[i] - y_ref[i];
            }
            diff = Rlange(&normtype, (mplapackint)n, (mplapackint)1, z, 1, dummywork);
            diffr = cast2double(diff);
            printf("         n       MFLOPS      error\n");
            printf("%10d   %10.3f   %10.3f\n", (int)n, (2.0 * (double)n) / elapsedtime * MFLOPS, diffr);
        } else {
            for (i = 0; i < n; i++) {
                x[i] = randomnumber(dummy);
                y[i] = y_ref[i] = randomnumber(dummy);
            }
            alpha = randomnumber(dummy);
            elapsedtime = 0.0;
            for (int j = 0; j < LOOPS; j++) {
                auto t1 = Clock::now();
                Raxpy(n, alpha, x, incx, y, incy);
                auto t2 = Clock::now();
                elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            }
            elapsedtime = elapsedtime / (double)LOOPS;
            printf("         n       MFLOPS     LOOPS\n");
            printf("%10d   %10.3f        %d\n", (int)n, (2.0 * (double)n) / elapsedtime * MFLOPS, (int)LOOPS);
        }
        delete[] z;
        delete[] y_ref;
        delete[] y;
        delete[] x;
        n = n + STEP;
    }
    if (check_flag)
        dlclose(handle);
}
