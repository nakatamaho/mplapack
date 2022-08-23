/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemv_qd.cpp,v 1.3 2010/08/07 05:50:09 nakatamaho Exp $
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
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_benchmark.h>

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_gemv(mplapackint m_i, mplapackint n_i) {
    double adds, muls, flops;
    double n, m;
    n = (double)n_i;
    m = (double)m_i;
    muls = m * n + 2. * m;
    adds = m * n;
    flops = muls + adds;
    return flops;
}

int main(int argc, char *argv[]) {
    mplapackint k, l, m, n;
    mplapackint STEPN = 7, STEPM = 7, N0, M0, LOOP = 3, TOTALSTEPS = 283;
    REAL alpha, beta, dummy, *dummywork;
    REAL mOne = -1;
    double elapsedtime, t1, t2;
    int i, p;
    int check_flag = 1;
    char trans = 'n', normtype = 'm';

    ___MPLAPACK_INITIALIZE___

    const char mpblas_sym[] = SYMBOL_GCC_RGEMV;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mpblas_ref)(const char *, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint);
    void (*raxpy_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    // initialization
    N0 = M0 = 1;
    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-M", argv[i]) == 0) {
                M0 = atoi(argv[++i]);
            } else if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STEPM", argv[i]) == 0) {
                STEPM = atoi(argv[++i]);
            } else if (strcmp("-T", argv[i]) == 0) {
                trans = 't';
            } else if (strcmp("-NOCHECK", argv[i]) == 0) {
                check_flag = 0;
            } else if (strcmp("-LOOP", argv[i]) == 0) {
                LOOP = atoi(argv[++i]);
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
        mpblas_ref = (void (*)(const char *, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint))dlsym(handle, mpblas_sym);
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

    n = N0;
    m = M0;
    for (p = 0; p < TOTALSTEPS; p++) {
        if (Mlsame(&trans, "n")) {
            k = n;
            l = m;
        } else {
            k = m;
            l = n;
        }
        REAL *x = new REAL[k];
        REAL *y = new REAL[l];
        REAL *yref = new REAL[l];
        REAL *A = new REAL[n * m];
        if (check_flag) {
            for (i = 0; i < k; i++) {
                x[i] = randomnumber(dummy);
            }
            for (i = 0; i < l; i++) {
                y[i] = yref[i] = randomnumber(dummy);
            }
            for (i = 0; i < k * l; i++) {
                A[i] = randomnumber(dummy);
            }
            alpha = randomnumber(dummy);
            beta = randomnumber(dummy);
            t1 = gettime();
            Rgemv(&trans, m, n, alpha, A, m, x, (mplapackint)1, beta, y, (mplapackint)1);
            t2 = gettime();
            elapsedtime = (t2 - t1);
            (*mpblas_ref)(&trans, m, n, alpha, A, m, x, (mplapackint)1, beta, yref, (mplapackint)1);
            (*raxpy_ref)(l, mOne, y, (mplapackint)1, yref, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)l, (mplapackint)1, yref, 1, dummywork);
            diffr = cast2double(diff);
            printf("     m       n      MFLOPS      error    trans\n");
            printf("%6d  %6d  %10.3f   %5.2e        %c\n", (int)n, (int)m, (2.0 * (double)n * (double)m) / elapsedtime * MFLOPS, diffr, trans);
        } else {
            for (i = 0; i < k; i++) {
                x[i] = randomnumber(dummy);
            }
            for (i = 0; i < l; i++) {
                y[i] = yref[i] = randomnumber(dummy);
            }
            for (i = 0; i < k * l; i++) {
                A[i] = randomnumber(dummy);
            }
            alpha = randomnumber(dummy);
            beta = randomnumber(dummy);
            elapsedtime = 0.0;
            for (int j = 0; j < LOOP; j++) {
                t1 = gettime();
                Rgemv(&trans, m, n, alpha, A, m, x, (mplapackint)1, beta, y, (mplapackint)1);
                t2 = gettime();
                elapsedtime = elapsedtime + (t2 - t1);
            }
            elapsedtime = elapsedtime / (double)LOOP;
            printf("     m       n      MFLOPS  trans\n");
            printf("%6d  %6d  %10.3f      %c\n", (int)n, (int)m, flops_gemv(m, n) / elapsedtime * MFLOPS, trans);
        }
        delete[] yref;
        delete[] y;
        delete[] x;
        delete[] A;
        n = n + STEPN;
        m = m + STEPM;
    }
    if (check_flag)
        dlclose(handle);
}
