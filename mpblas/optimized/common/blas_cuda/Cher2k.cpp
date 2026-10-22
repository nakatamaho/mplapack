/*
 * Copyright (c) 2026
 *      Nakata, Maho
 *      All rights reserved.
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

// Cher2k of the dd/td/qd CUDA libraries: eligible calls run on the GPU
// (Cblas_bridge_cuda.cpp, ../Cblas_cuda_common.cuh), the others on the CPU
// (../blas/Cher2k_omp.cpp, then the reference code).  All of them give the
// reference result bit for bit.

#include <mpblas.h>

bool Cher2k_gpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc);
bool Cher2k_omp(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc);
void Cher2k_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc);

void Cher2k(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc)
{
    if (Cher2k_gpu(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc))
        return;
    if (Cher2k_omp(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc))
        return;
    Cher2k_reference_cpu(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
