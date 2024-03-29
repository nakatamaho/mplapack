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

#include <iostream>
#include <stdio.h>
#include "dd_real_cuda.h"
#include <mplapack_config.h>
#include <cuda_runtime.h>
#include <cuda.h>

mplapackint Mlsame_dd(const char *a, const char *b);
void Mxerbla_dd(const char *srname, int info);

__global__ void Rsyrk_NU_0 (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_NU_p (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TU_0 (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TU_p (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_NL_0 (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_NL_p (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TL_0 (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);
__global__ void Rsyrk_TL_p (dd_real * Adev, dd_real * Cdev, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc, dd_real alpha, dd_real beta);

// matrix block size
#define Bk  (16)
#define Bn  (16)
#define Gn   (4)

extern void Is_cuda_Rgemm_error(cudaError_t rc, const char *mes, mplapackint n, mplapackint k, mplapackint lda, mplapackint ldc);

void Rsyrk_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, dd_real alpha, dd_real * Adev, mplapackint lda, dd_real beta, dd_real * Cdev, mplapackint ldc)
{
    mplapackint upper;
    cudaError_t rc;
    dd_real Zero, One;

    dd_set(Zero, 0.0, 0.0);
    dd_set(One, 1.0, 0.0);

    //quick return if possible.
    if ((n == 0) || (( dd_eq(alpha, Zero) || (k == 0)) && dd_eq(beta, One) ))
	return;

    upper = Mlsame_dd(uplo, "U");

    //start the operations.
    if (Mlsame_dd(trans, "N")) {
	//Form C := alpha*A*A' + beta*C.
	if (upper) {
            // calculating and updating C
            dim3 grid(n / Bn + (n % Bn == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bn, Bn);
            if(k % Bk == 0 && n % (Gn * Bn) == 0){
                Rsyrk_NU_0 <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }else{
                Rsyrk_NU_p <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }
	} else {
            // calculating and updating C
            dim3 grid(n / Bn + (n % Bn == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bn, Bn);
            if(k % Bk == 0 && n % (Gn * Bn) == 0){
                Rsyrk_NL_0 <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }else{
                Rsyrk_NL_p <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }
	}
    } else {
	//Form C := alpha*A'*A + beta*C.
	if (upper) {
            // calculating and updating C
            dim3 grid(n / Bn + (n % Bn == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bn, Bn);
            if(k % Bk == 0 && n % (Gn * Bn) == 0){
                Rsyrk_TU_0 <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }else{
                Rsyrk_TU_p <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }
	} else {
            // calculating and updating C
            dim3 grid(n / Bn + (n % Bn == 0 ? 0 : 1), n / (Gn * Bn)
                      + (n % (Gn * Bn) == 0 ? 0 : 1)), block(Bn, Bn);
            if(k % Bk == 0 && n % (Gn * Bn) == 0){
                Rsyrk_TL_0 <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }else{
                Rsyrk_TL_p <<< grid, block >>> (Adev, Cdev, n, k, lda, ldc, alpha, beta);
            }
	}
    }
    cudaThreadSynchronize();
}

