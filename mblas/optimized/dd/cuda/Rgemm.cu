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

void Rgemm_fermi(const char *transa, const char *transb, mpackint m, mpackint n, mpackint k, dd_real alpha, dd_real * A, mpackint lda, dd_real * B, mpackint ldb, dd_real beta, dd_real * C, mpackint ldc);

void Rgemm_tesla(const char *transa, const char *transb, mpackint m, mpackint n, mpackint k, dd_real alpha, dd_real * A, mpackint lda, dd_real * B, mpackint ldb, dd_real beta, dd_real * C, mpackint ldc);

void Rgemm(const char *transa, const char *transb, mpackint m, mpackint n, mpackint k, dd_real alpha, dd_real * A, mpackint lda, dd_real * B, mpackint ldb, dd_real beta, dd_real * C, mpackint ldc)
{
//just call fermi version
    Rgemm_fermi(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

//if tesla version
//  Rgemm_tesla(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    return;
}
