/*
 * Copyright (c) 2010-2011
 *      RIKEN
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
/*
  Contributed by Takao, Yasuyoshi and Nakata, Maho, 2010-2011
*/
/*
Based on http://www.netlib.org/blas/dgemm.f
Rgemm performs one of the matrix-matrix operations
 C := alpha*op(A)*op(B) + beta*C,
where op(X) is one of
 op(X) = X or op(X) = X',
alpha and beta are scalars, and A, B and C are matrices, with op( A )
an m by k matrix, op(B) a k by n matrix and C an m by n matrix.
*/

#include <iostream>
#include <stdio.h>
#include "dd_real_cuda.h"
#include <mpack_config.h>
#include <cuda_runtime.h>
#include <cuda.h>

mpackint Mlsame_dd(const char *a, const char *b);
void Mxerbla_dd(const char *srname, int info);

// define texture memory
texture < int4, 1 > tex_x_double_A;
texture < int4, 1 > tex_x_double_B;

// matrix block size
#define Bk  (16)
#define Bn  (16)
#define Gn   (4)

static
__inline__ __device__ dd_real fetch_x_A(const int &i)
{
    register int4 v = tex1Dfetch(tex_x_double_A, i);
    register
    dd_real r;
    r.x[0] = __hiloint2double(v.y, v.x);
    r.x[1] = __hiloint2double(v.w, v.z);
    return r;
}

__global__ void Rsyrk_NU_0 (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_NU_p (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TU_0 (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TU_p (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_NL_0 (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_NL_p (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TL_0 (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TL_p (dd_real * Adev, dd_real * Cdev, mpackint n, mpackint k, mpackint lda, mpackint ldc, dd_real alpha, dd_real beta);

void Rsyrk_cuda(const char *uplo, const char *trans, mpackint n, mpackint k, dd_real alpha, dd_real * Adev, mpackint lda, dd_real beta, dd_real * Cdev, mpackint ldc);

void Is_cuda_Rgemm_error(cudaError_t rc, const char *mes, mpackint n, mpackint k, mpackint lda, mpackint ldc)
{
    if (rc != cudaSuccess) {
	fprintf(stderr, "%s : n:%d k:%d lda:%d ldc:%d\n", mes, n, k, lda, ldc);
	exit(1);
    }
    /* not an error */
}

#include <Rsyrk_NU_0.cu>
#include <Rsyrk_NU_p.cu>
#include <Rsyrk_TU_0.cu>
#include <Rsyrk_TU_p.cu>
#include <Rsyrk_NL_0.cu>
#include <Rsyrk_NL_p.cu>
#include <Rsyrk_TL_0.cu>
#include <Rsyrk_TL_p.cu>

void Rsyrk(const char *uplo, const char *trans, mpackint n, mpackint k, dd_real alpha, dd_real * A, mpackint lda, dd_real beta, dd_real * C, mpackint ldc)
{
    mpackint i, j, nota, upper, info, nrowa;
    dd_real temp, Zero, One;
    cudaError_t rc;

    dd_set(Zero, 0.0, 0.0);
    dd_set(One, 1.0, 0.0);

    //Test the input parameters.
    if (Mlsame_dd(trans, "N"))
	nrowa = n;
    else
	nrowa = k;
    nota = Mlsame_dd(trans, "N");
    upper = Mlsame_dd(uplo, "U");

    info = 0;
    if ((!upper) && (!Mlsame_dd(uplo, "L")))
	info = 1;
    else if ((!Mlsame_dd(trans, "N")) && (!Mlsame_dd(trans, "T"))
	&& (!Mlsame_dd(trans, "C")))
	info = 2;
    else if (n < 0)
	info = 3;
    else if (k < 0)
	info = 4;
    else if (lda < std::max((mpackint) 1, nrowa))
	info = 7;
    else if (ldc < std::max((mpackint) 1, n))
	info = 10;
    if (info != 0) {
	Mxerbla_dd("Rsyrk ", info);
	return;
    }

    //quick return if possible.
    if ((n == 0) || ((dd_eq(alpha, Zero) || (k == 0)) && dd_eq(beta, One)))
	return;

    //allocate device memory for GPU
    dd_real *Adev, *Cdev;
    int size_A, size_C;
    if (nota)
	size_A = lda * k - (lda - n);
    else
	size_A = lda * n - (lda - k);
    size_C = ldc * n - (ldc - n);
    rc = cudaMalloc((void **) &Adev, size_A * sizeof(dd_real));
        Is_cuda_Rgemm_error(rc, "cudaMalloc A error", n, k, lda, ldc);
    rc = cudaMalloc((void **) &Cdev, size_C * sizeof(dd_real));
        Is_cuda_Rgemm_error(rc, "cudaMalloc C error", n, k, lda, ldc);
    rc = cudaMemcpy(Adev, A, size_A * sizeof(dd_real), cudaMemcpyHostToDevice);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy A error", n, k, lda, ldc);
    rc = cudaMemcpy(Cdev, C, size_C * sizeof(dd_real), cudaMemcpyHostToDevice);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy C error", n, k, lda, ldc);

//And when alpha == 0.0
    if (dd_eq(alpha, Zero)) {
	if (dd_eq(beta, Zero)) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
		    C[i + j * ldc] = Zero;
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
		    dd_mul_host(beta, C[i + j * ldc], C[i + j * ldc]);
		}
	    }
	}
	return;
    }

    Rsyrk_cuda(uplo, trans, n, k, alpha, Adev, lda, beta, Cdev, ldc);

    rc = cudaMemcpy(C, Cdev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy C error", n, k, lda, ldc);
    rc = cudaFree(Adev);
        Is_cuda_Rgemm_error(rc, "cudaFree A error", n, k, lda, ldc);
    rc = cudaFree(Cdev);
        Is_cuda_Rgemm_error(rc, "cudaFree C error", n, k, lda, ldc);
    return;
}
