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

#include <cuda.h>
#include <cuda_runtime.h>

void Rgemm_cuda(const char *transa, const char *transb, mplapackint m, mplapackint n, mplapackint k, REAL alpha, REAL *Adev, mplapackint lda, REAL *Bdev, mplapackint ldb, REAL beta, REAL *Cdev, mplapackint ldc);

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

void SetDevice() {
    int gpudevice = 0; // cpu number 0,1,...

    int device_count = 0;
    int device;

    cudaGetDeviceCount(&device_count);
    printf("device_count : %d\n", device_count);

    cudaError_t cudareturn;
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, gpudevice);

    printf("device name -> %s \n", deviceProp.name);

    if (deviceProp.warpSize <= 1) {
        printf("warning, CUDA Device Emulation (CPU)  detected, exiting\n");
        exit(1);
    }
    // set GPU that is to be used
    cudareturn = cudaSetDevice(gpudevice);
    printf("cudareturn -> %d\n", cudareturn);

    if (cudareturn == cudaErrorInvalidDevice) {
        perror("cudaSetDevice returned  cudaErrorInvalidDevice");
    } else {
        cudaGetDevice(&device);
        printf("cudaGetDevice()=%d\n", device);
    }
}

int main(int argc, char *argv[]) {
    REAL alpha, beta, dummy;
    REAL *dummywork;
    double elapsedtime;
    double *dummyd;
    char transa, transb, normtype;
    int N0, M0, K0, STEPN = 7, STEPM = 7, STEPK = 7, LOOP = 3, TOTALSTEPS = 720;
    int lda, ldb, ldc;
    int i, j, m, n, k, ka, kb, p, q;
    int check_flag = 1;

    const char mpblas_sym[] = SYMBOL_GCC_RGEMM;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mpblas_ref)(const char *, const char *, mplapackint, mplapackint, mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint, REAL, REAL *, mplapackint);
    void (*raxpy_ref)(mplapackint, REAL, REAL *, mplapackint, REAL *, mplapackint);
    char *error;
    REAL diff;
    double diffr;

    using Clock = std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;

    // initialization
    N0 = M0 = K0 = 1;
    STEPM = STEPN = STEPK = 1;
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

    SetDevice();
    // dummy memory allocation for initialization
    cudaMalloc((void **)&dummyd, 16);
    cudaFree(dummyd);

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

        REAL *A = new REAL[lda * ka];
        REAL *B = new REAL[ldb * kb];
        REAL *C = new REAL[ldc * n];
        REAL *Cref = new REAL[ldc * n];
        REAL mOne = -1;
        alpha = randomnumber(dummy);
        beta = randomnumber(dummy);
        for (i = 0; i < lda * ka; i++) {
            A[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldb * kb; i++) {
            B[i] = randomnumber(dummy);
        }
        for (i = 0; i < ldc * n; i++) {
            C[i] = Cref[i] = randomnumber(dummy);
        }
        dd_real *Adev, *Bdev, *Crefev;
        int size_A, size_B, size_C;
        if (Mlsame(&transa, "n")) {
            size_A = lda * k - (lda - m);
        } else
            size_A = lda * m - (lda - k);
        if (Mlsame(&transb, "n")) {
            size_B = ldb * n - (ldb - k);
        } else
            size_B = ldb * k - (ldb - n);
        size_C = ldc * n - (ldc - m);

        if (check_flag) {
            cudaMalloc((void **)&Adev, size_A * sizeof(REAL));
            cudaMalloc((void **)&Bdev, size_B * sizeof(REAL));
            cudaMalloc((void **)&Crefev, size_C * sizeof(REAL));
            cudaMemcpy(Adev, A, size_A * sizeof(REAL), cudaMemcpyHostToDevice);
            cudaMemcpy(Bdev, B, size_B * sizeof(REAL), cudaMemcpyHostToDevice);
            cudaMemcpy(Crefev, C, size_C * sizeof(REAL), cudaMemcpyHostToDevice);

            auto t1 = Clock::now();
            Rgemm_cuda(&transa, &transb, m, n, k, alpha, Adev, lda, Bdev, ldb, beta, Crefev, ldc);
            auto t2 = Clock::now();
            cudaMemcpy(C, Crefev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
            cudaFree(Adev);
            cudaFree(Bdev);
            cudaFree(Crefev);
            elapsedtime = (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;

            (*mpblas_ref)(&transa, &transb, m, n, k, alpha, A, lda, B, ldb, beta, Cref, ldc);
            (*raxpy_ref)((mplapackint)(ldc * n), mOne, C, (mplapackint)1, Cref, (mplapackint)1);
            diff = Rlange(&normtype, (mplapackint)ldc, (mplapackint)n, Cref, ldc, dummywork);
            diffr = cast2double(diff);
            printf("    m     n     k     MFLOPS      error   transa   transb\n");
            printf("%5d %5d %5d  %10.3f    %5.2e     %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, diffr, transa, transb);
        } else {
            cudaMalloc((void **)&Adev, size_A * sizeof(REAL));
            cudaMalloc((void **)&Bdev, size_B * sizeof(REAL));
            cudaMalloc((void **)&Crefev, size_C * sizeof(REAL));
            cudaMemcpy(Adev, A, size_A * sizeof(REAL), cudaMemcpyHostToDevice);
            cudaMemcpy(Bdev, B, size_B * sizeof(REAL), cudaMemcpyHostToDevice);
            cudaMemcpy(Crefev, C, size_C * sizeof(REAL), cudaMemcpyHostToDevice);

            elapsedtime = 0.0;
            for (int j = 0; j < LOOP; j++) {
                auto t1 = Clock::now();
                Rgemm_cuda(&transa, &transb, m, n, k, alpha, Adev, lda, Bdev, ldb, beta, Crefev, ldc);
                auto t2 = Clock::now();
                elapsedtime = elapsedtime + (double)duration_cast<nanoseconds>(t2 - t1).count() / 1.0e9;
            }
            elapsedtime = elapsedtime / (double)LOOP;

            cudaMemcpy(C, Crefev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
            cudaFree(Adev);
            cudaFree(Bdev);
            cudaFree(Crefev);
            printf("    m     n     k     MFLOPS    transa   transb\n");
            printf("%5d %5d %5d  %10.3f        %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, transa, transb);
        }
        delete[] Cref;
        delete[] C;
        delete[] B;
        delete[] A;
        m = m + STEPM;
        n = n + STEPN;
        k = k + STEPK;
    }
    if (check_flag)
        dlclose(handle);
}
