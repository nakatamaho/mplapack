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

// Rtrmm of the dd/td/qd CUDA libraries: eligible calls run on the GPU
// (Rblas_bridge_cuda.cpp, ../Rblas_cuda_common.cuh), the others on the CPU
// (../blas/Rtrmm_omp.cpp: lane-parallel OpenMP, else the reference code).
// All of them give the reference result bit for bit.

#include <mpblas.h>

bool Rtrmm_gpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb);
bool Rtrmm_omp(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb);
void Rtrmm_reference_cpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb);

void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb)
{
    if (Rtrmm_gpu(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb))
        return;
    if (Rtrmm_omp(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb))
        return;
    Rtrmm_reference_cpu(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}
