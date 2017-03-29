/*
 * Copyright (c) 2008-2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemm_dd.cpp,v 1.4 2010/08/07 05:50:09 nakatamaho Exp $
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

#define ___MPACK_BUILD_WITH_DD___
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include <mblas.h>
#include <mlapack.h>
#include <mpack_benchmark.h>

#include <cuda.h>
#include <cuda_runtime.h>

#define TOTALSTEPS 1000

void Rsyrk_cuda(const char *uplo, const char *trans, mpackint n, mpackint k, dd_real alpha, dd_real * Adev, mpackint lda, dd_real beta, dd_real * Cdev, mpackint ldc);

void SetDevice()
{
    int gpudevice = 0;		// cpu number 0,1,...
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
    //set GPU that is to be used
    cudareturn = cudaSetDevice(gpudevice);
    printf("cudareturn -> %d\n", cudareturn);

    if (cudareturn == cudaErrorInvalidDevice) {
	perror("cudaSetDevice returned  cudaErrorInvalidDevice");
    } else {
	cudaGetDevice(&device);
	printf("cudaGetDevice()=%d\n", device);
    }
}

int main(int argc, char *argv[])
{
    REAL alpha, beta, dummy;
    REAL *dummywork;
    double elapsedtime, t1, t2;
    double *dummyd;
    char uplo, trans, normtype;
    int N0, K0, STEPN, STEPK;
    int lda, ldc;
    int i, j, n, k, ka, kb, p, q;
    int check_flag = 1;

    const char mblas_sym[] = SYMBOL_GCC_RSYRK;
    const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
    void *handle;
    void (*mblas_ref) (const char *, const char *, mpackint, mpackint, REAL, REAL *, mpackint, REAL, REAL *, mpackint);
    void (*raxpy_ref) (mpackint, REAL, REAL *, mpackint, REAL *, mpackint);
    char *error;
    REAL diff;
    double diffr;

    // initialization
    N0 = K0 = 1;
    STEPN = STEPK = 1;
    uplo = 'u';
    trans = 'n';
    normtype = 'm';
    if (argc != 1) {
	for (i = 1; i < argc; i++) {
	    if (strcmp("-N", argv[i]) == 0) {
		N0 = atoi(argv[++i]);
	    } else if (strcmp("-K", argv[i]) == 0) {
		K0 = atoi(argv[++i]);
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
	    }
	}
    }
    if (check_flag) {
	handle = dlopen(MBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
	if (!handle) {
	    printf("dlopen: %s\n", dlerror());
	    return 1;
	}
	mblas_ref = (void (*)
		     (const char *, const char *, mpackint, mpackint, REAL, REAL *, mpackint, REAL, REAL *, mpackint)) dlsym(handle, mblas_sym);
	if ((error = dlerror()) != NULL) {
	    fprintf(stderr, "%s\n", error);
	    return 1;
	}

	raxpy_ref = (void (*)(mpackint, REAL, REAL *, mpackint, REAL *, mpackint))
	    dlsym(handle, raxpy_sym);
	if ((error = dlerror()) != NULL) {
	    fprintf(stderr, "%s\n", error);
	    return 1;
	}
    }

    SetDevice();
    //dummy memory allocation for initialization
    cudaMalloc((void **) &dummyd, 16);
    cudaFree(dummyd);

    n = N0;
    k = K0;
    for (p = 0; p < TOTALSTEPS; p++) {
	if (Mlsame(&trans, "n")) {
	    ka = k;
	    lda = n;
	} else {
	    ka = n;
	    lda = k;
	}
	ldc = n;

	REAL *A = new REAL[lda * ka];
	REAL *C = new REAL[ldc * n];
	REAL *Cd = new REAL[ldc * n];
	REAL mOne = -1;
	alpha = randomnumber(dummy);
	beta = randomnumber(dummy);
	for (i = 0; i < lda * ka; i++) {
	    A[i] = randomnumber(dummy);
	}
	for (i = 0; i < ldc * n; i++) {
	    C[i] = Cd[i] = randomnumber(dummy);
	}
	dd_real *Adev, *Cdev;
	int size_A, size_C;
	if (Mlsame(&trans, "n")) {
	    size_A = lda * k - (lda - n);
	} else
	    size_A = lda * n - (lda - k);
	size_C = ldc * n - (ldc - n);

	if (check_flag) {
	    cudaMalloc((void **) &Adev, size_A * sizeof(REAL));
	    cudaMalloc((void **) &Cdev, size_C * sizeof(REAL));
	    cudaMemcpy(Adev, A, size_A * sizeof(REAL), cudaMemcpyHostToDevice);
	    cudaMemcpy(Cdev, C, size_C * sizeof(REAL), cudaMemcpyHostToDevice);

	    t1 = gettime();
	    Rsyrk_cuda(&uplo, &trans, n, k, alpha, Adev, lda, beta, Cdev, ldc);
	    t2 = gettime();
	    cudaMemcpy(C, Cdev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
	    cudaFree(Adev);
	    cudaFree(Cdev);

	    elapsedtime = (t2 - t1);
	    (*mblas_ref) (&uplo, &trans, n, k, alpha, A, lda, beta, Cd, ldc);
	    (*raxpy_ref) ((mpackint) (ldc * n), mOne, C, (mpackint) 1, Cd, (mpackint) 1);
	    diff = Rlange(&normtype, (mpackint) ldc, (mpackint) n, Cd, ldc, dummywork);
	    diffr = cast2double(diff);
            printf ("    n     k      MFLOPS       error    uplo    trans\n");
	    //2n^2k+2n^2 flops are needed
	    printf("%5d %5d  %10.3f    %5.2e       %c        %c\n", (int) n, (int) k, (2.0 * (double) n * (double) n * (double) k + 2.0 * (double) n * (double) n) / elapsedtime * MFLOPS, diffr, uplo, trans);
	} else {
	    cudaMalloc((void **) &Adev, size_A * sizeof(REAL));
	    cudaMalloc((void **) &Cdev, size_C * sizeof(REAL));
	    cudaMemcpy(Adev, A, size_A * sizeof(REAL), cudaMemcpyHostToDevice);
	    cudaMemcpy(Cdev, C, size_C * sizeof(REAL), cudaMemcpyHostToDevice);

	    t1 = gettime();
	    Rsyrk_cuda(&uplo, &trans, n, k, alpha, Adev, lda, beta, Cdev, ldc);
	    t2 = gettime();

	    cudaMemcpy(C, Cdev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
	    cudaFree(Adev);
	    cudaFree(Cdev);

	    elapsedtime = (t2 - t1);
            printf ("    n     k      MFLOPS       error    uplo    trans\n");
	    //2n^2k+2n^2 flops are needed
	    printf("%5d %5d %10.3f     %5.2e        %c\n", (int) n, (int) k, (2.0 * (double) n * (double) n * (double) k + 2.0 * (double) n * (double) n) / elapsedtime * MFLOPS, uplo, trans);
	}
	delete[]Cd;
	delete[]C;
	delete[]A;
	n = n + STEPN;
	k = k + STEPK;
    }
    if (check_flag)
	dlclose(handle);
}
