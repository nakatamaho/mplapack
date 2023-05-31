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
#include <mplapack_config.h>
#include <cuda_runtime.h>
#include <cuda.h>

mplapackint Mlsame_dd(const char *a, const char *b);
void Mxerbla_dd(const char *srname, int info);

// matrix block size
#define Bm  (16)
#define Bk  (16)
#define Bn  (16)
#define Gn   (2)

//for alpha*A*B + beta
__global__ void Rgemm_tesla_NN_0 (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rgemm_tesla_NN_p (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);

//for alpha*A'*B + beta
__global__ void Rgemm_tesla_TN_0 (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rgemm_tesla_TN_p (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);

//for alpha*A*B' + beta
__global__ void Rgemm_tesla_NT_0 (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rgemm_tesla_NT_p (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);

//for alpha*A'*B' + beta
__global__ void Rgemm_tesla_TT_0 (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rgemm_tesla_TT_p (dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc, dd_real alpha, dd_real beta);

void Is_cuda_Rgemm_error(cudaError_t rc, const char *mes, mplapackint m, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldb, mplapackint ldc);

#define fetch_x_A(XX,i) XX[i]
#define fetch_x_B(XX,i) XX[i]

#include <Rgemm_tesla_NN_0.cu>
#include <Rgemm_tesla_NN_p.cu>
#include <Rgemm_tesla_TN_0.cu>
#include <Rgemm_tesla_TN_p.cu>
#include <Rgemm_tesla_NT_0.cu>
#include <Rgemm_tesla_NT_p.cu>
#include <Rgemm_tesla_TT_0.cu>
#include <Rgemm_tesla_TT_p.cu>

void Rgemm_tesla_cuda(const char *transa, const char *transb, mplapackint m, mplapackint n, mplapackint k, dd_real alpha, dd_real * Adev, mplapackint lda, dd_real * Bdev, mplapackint ldb, dd_real beta, dd_real * Cdev, mplapackint ldc)
{
    mplapackint nota, notb;
    cudaError_t rc;

    nota = Mlsame_dd(transa, "N");
    notb = Mlsame_dd(transb, "N");

    if (notb) {
        if (nota) {
	    //Form C := alpha*A*B + beta*C.
            // calculating and updating C
            dim3 grid(m / Bm + (m % Bm == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bm, Bn);
            if(m % Bm == 0 && k % Bk == 0 && n % (Gn * Bn) == 0){
                Rgemm_tesla_NN_0 <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }else{
                Rgemm_tesla_NN_p <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }
        } else {
	    //Form C := alpha*A'*B + beta*C.
            // calculating and updating C
            dim3 grid(m / Bm + (m % Bm == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bm, Bn);
            if(m % Bm == 0 && k % Bk == 0 && n % (Gn * Bn) == 0){
                Rgemm_tesla_TN_0 <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }else{
                Rgemm_tesla_TN_p <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }
        }
    } else {
        if (nota) {
	    //Form  C := alpha*A*B' + beta*C.
            // calculating and updating C
            dim3 grid(m / Bm + (m % Bm == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bm, Bn);
            if(m % Bm == 0 && k % Bk == 0 && n % (Gn * Bn) == 0){
                Rgemm_tesla_NT_0 <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }else{
                Rgemm_tesla_NT_p <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }
        } else {
	    //Form  C := alpha*A'*B' + beta*C.
            // calculating and updating C
            dim3 grid(m / Bm + (m % Bm == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bm, Bn);
            if(m % Bm == 0 && k % Bk == 0 && n % (Gn * Bn) == 0){
                Rgemm_tesla_TT_0 <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }else{
                Rgemm_tesla_TT_p <<< grid, block >>> (Adev, Bdev, Cdev, m, n, k, lda, ldb, ldc, alpha, beta);
            }
        }
    }
    cudaThreadSynchronize();
}

void Rgemm_tesla(const char *transa, const char *transb, mplapackint m, mplapackint n, mplapackint k, dd_real alpha, dd_real * A, mplapackint lda, dd_real * B, mplapackint ldb, dd_real beta, dd_real * C, mplapackint ldc)
{
    mplapackint i, j, nota, notb, nrowa, nrowb, ncola, info;
    dd_real temp, Zero, One;
    cudaError_t rc;

    dd_set(Zero, 0.0, 0.0);
    dd_set(One, 1.0, 0.0);

    nota = Mlsame_dd(transa, "N");
    notb = Mlsame_dd(transb, "N");
    if (nota) {
	nrowa = m;
	ncola = k;
    } else {
	nrowa = k;
	ncola = m;
    }
    if (notb) {
	nrowb = k;
    } else {
	nrowb = n;
    }
//Test the input parameters.
    info = 0;
    if (!nota && (!Mlsame_dd(transa, "C"))
	&& (!Mlsame_dd(transa, "T")))
	info = 1;
    else if (!notb && (!Mlsame_dd(transb, "C"))
	     && (!Mlsame_dd(transb, "T")))
	info = 2;
    else if (m < 0)
	info = 3;
    else if (n < 0)
	info = 4;
    else if (k < 0)
	info = 5;
    else if (lda < std::max((mplapackint) 1, nrowa))
	info = 8;
    else if (ldb < std::max((mplapackint) 1, nrowb))
	info = 10;
    else if (ldc < std::max((mplapackint) 1, m))
	info = 13;
    if (info != 0) {
	Mxerbla_dd("Rgemm ", info);
	return;
    }
    //Quick return if possible.
    if ((m == 0)
	|| (n == 0)
	|| ((dd_eq(alpha, Zero)
	     || (k == 0))
	    && dd_eq(beta, One))) {
	return;
    }

    //allocate device memory for GPU
    dd_real *Adev, *Bdev, *Cdev;
    int size_A, size_B, size_C;
    if (nota)
	size_A = lda * k - (lda - m);
    else
	size_A = lda * m - (lda - k);
    if (notb)
	size_B = ldb * n - (ldb - k);
    else
	size_B = ldb * k - (ldb - n);
    size_C = ldc * n - (ldc - m);
    rc = cudaMalloc((void **) &Adev, size_A * sizeof(dd_real));
        Is_cuda_Rgemm_error(rc, "cudaMalloc A error", m, n, k, lda, ldb, ldc);
    rc = cudaMalloc((void **) &Bdev, size_B * sizeof(dd_real));
        Is_cuda_Rgemm_error(rc, "cudaMalloc B error", m, n, k, lda, ldb, ldc);
    rc = cudaMalloc((void **) &Cdev, size_C * sizeof(dd_real));
        Is_cuda_Rgemm_error(rc, "cudaMalloc C error", m, n, k, lda, ldb, ldc);
    rc = cudaMemcpy(Adev, A, size_A * sizeof(dd_real), cudaMemcpyHostToDevice);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy A error", m, n, k, lda, ldb, ldc);
    rc = cudaMemcpy(Bdev, B, size_B * sizeof(dd_real), cudaMemcpyHostToDevice);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy B error", m, n, k, lda, ldb, ldc);
    rc = cudaMemcpy(Cdev, C, size_C * sizeof(dd_real), cudaMemcpyHostToDevice);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy C error", m, n, k, lda, ldb, ldc);

//And when alpha == 0.0
    if (dd_eq(alpha, Zero)) {
	if (dd_eq(beta, Zero)) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = Zero;
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    dd_mul_host(beta, C[i + j * ldc], C[i + j * ldc]);
		}
	    }
	}
	return;
    }

    Rgemm_tesla_cuda(transa, transb, m, n, k, alpha, Adev, lda, Bdev, ldb, beta, Cdev, ldc);

    rc = cudaMemcpy(C, Cdev, size_C * sizeof(dd_real), cudaMemcpyDeviceToHost);
        Is_cuda_Rgemm_error(rc, "cudaMemcpy C error", m, n, k, lda, ldb, ldc);
    rc = cudaFree(Adev);
        Is_cuda_Rgemm_error(rc, "cudaFree A error", m, n, k, lda, ldb, ldc);
    rc = cudaFree(Bdev);
        Is_cuda_Rgemm_error(rc, "cudaFree B error", m, n, k, lda, ldb, ldc);
    rc = cudaFree(Cdev);
        Is_cuda_Rgemm_error(rc, "cudaFree C error", m, n, k, lda, ldb, ldc);
    return;
}
