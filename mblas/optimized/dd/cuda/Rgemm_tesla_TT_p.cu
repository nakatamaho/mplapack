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

//for alpha*A'*B' + beta*C
__global__ void Rgemm_tesla_TT_p(dd_real * Adev, dd_real * Bdev, dd_real * Cdev, mpackint m, mpackint n, mpackint k, mpackint lda, mpackint ldb, mpackint ldc, dd_real alpha, dd_real beta)
{
    int i, j;
    dd_real c_val0;
    dd_real c_val1;
    int iAb, jAb, A_i, A_j, Ab_i, Ab_j;
    int iBb, jBb, B_i, B_j, Bb_i, Bb_j;
    int iCb, jCb, C_i, C_j;
    dd_real Aval;
    int Lp;
    dd_real regA;
    dd_real regB0;
    dd_real regB1;
    dd_real temp0;
    dd_real temp1;

    __shared__ dd_real Ab[Bm][Bk + 1];
    __shared__ dd_real Bb[2][Bn][Bk + 1];

    Ab_j = threadIdx.x; //exchange x for y for coalescing
    Bb_i = threadIdx.y; //exchange x for y for coalescing

    iAb = blockIdx.x;
    A_i = blockDim.x * iAb + threadIdx.y; //exchange x for y for coalescing

    //load first data of A from global memory into register
    A_j = blockDim.y * 0 + threadIdx.x; //exchange x for y for coalescing
    regA = fetch_x_A(min(A_i, (int) (m - 1)) * lda + min(A_j, (int) (k - 1)));

    //load first data of B from global memory into register
    iBb = 0;
    B_i = blockDim.x * iBb + threadIdx.y; //exchange x for y for coalescing

    jBb = blockIdx.y * Gn + 0;
    B_j = blockDim.y * jBb + threadIdx.x; //exchange x for y for coalescing
    regB0 = fetch_x_B(min(B_i, (int) (k - 1)) * ldb + min(B_j, (int) (n - 1)));

    jBb = blockIdx.y * Gn + 1;
    B_j = blockDim.y * jBb + threadIdx.x; //exchange x for y for coalescing
    regB1 = fetch_x_B(min(B_i, (int) (k - 1)) * ldb + min(B_j, (int) (n - 1)));

    // get initial Cdev data
    iCb = blockIdx.x;
    C_i = blockDim.x * iCb + threadIdx.x;

    jCb = blockIdx.y * Gn + 0;
    C_j = blockDim.y * jCb + threadIdx.y;
    temp0 = Cdev[min(C_i, (int) (m - 1)) + min(C_j, (int) (n - 1)) * ldc];

    jCb = blockIdx.y * Gn + 1;
    C_j = blockDim.y * jCb + threadIdx.y;
    temp1 = Cdev[min(C_i, (int) (m - 1)) + min(C_j, (int) (n - 1)) * ldc];

    c_val0.x[0] = c_val0.x[1] = c_val1.x[0] = c_val1.x[1] = 0.0;

    for (i = 0; i < k / Bk + (k % Bk == 0 ? 0 : 1); i++) {
        Ab_i = threadIdx.y; //exchange x for y for coalescing
	// load data into Ab (in shared memory) from register
	Ab[Ab_i][Ab_j] = regA;

	// load data into Bb (in shared memory) from register
        Bb_j = threadIdx.x + 0 * blockDim.y; //exchange x for y for coalescing
	Bb[0][Bb_j][Bb_i] = regB0;
	Bb[1][Bb_j][Bb_i] = regB1;

	// syncronize in the block
	__syncthreads();

	// update C value
	Lp = (k % Bk == 0 ? Bk : (i < k / Bk ? Bk : k % Bk));

        Ab_i = threadIdx.x; //recover Ab_i for coalescing
        Bb_j = threadIdx.y + 0 * blockDim.y; //recover Bb_j for coalescing

#pragma unroll
	for (j = 0; j < Lp; j++) {
	    //take advantage of speed difference between register and smem
	    Aval = Ab[Ab_i][j];
	    dd_mad(c_val0, Aval, Bb[0][Bb_j][j]);
	    dd_mad(c_val1, Aval, Bb[1][Bb_j][j]);
	}

	// load next data of A from global memory into register
	jAb = i + 1;
        A_j = blockDim.y * jAb + threadIdx.x; //exchange x for y for coalescing
	regA = fetch_x_A(min(A_i, (int) (m - 1)) * lda + min(A_j, (int) (k - 1)));

	// load next data of B from global memory into register
	iBb = i + 1;
        B_i = blockDim.x * iBb + threadIdx.y; //exchange x for y for coalescing

	jBb = blockIdx.y * Gn + 0;
        B_j = blockDim.y * jBb + threadIdx.x; //exchange x for y for coalescing
	regB0 = fetch_x_B(min(B_i, (int) (k - 1)) * ldb + min(B_j, (int) (n - 1)));

	jBb = blockIdx.y * Gn + 1;
        B_j = blockDim.y * jBb + threadIdx.x; //exchange x for y for coalescing
	regB1 = fetch_x_B(min(B_i, (int) (k - 1)) * ldb + min(B_j, (int) (n - 1)));

	__syncthreads();
    }

    dd_mul(c_val0, alpha, c_val0);
    dd_mul(c_val1, alpha, c_val1);

    jCb = blockIdx.y * Gn + 0;
    C_j = blockDim.y * jCb + threadIdx.y;
    dd_mul(beta, temp0, temp0);
    dd_add(temp0, c_val0, Cdev[min(C_i, (int) (m - 1)) + min(C_j, (int) (n - 1)) * ldc]);

    jCb = blockIdx.y * Gn + 1;
    C_j = blockDim.y * jCb + threadIdx.y;
    dd_mul(beta, temp1, temp1);
    dd_add(temp1, c_val1, Cdev[min(C_i, (int) (m - 1)) + min(C_j, (int) (n - 1)) * ldc]);

}
