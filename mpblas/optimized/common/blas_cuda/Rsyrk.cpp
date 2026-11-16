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

// Rsyrk of the dd/td/qd CUDA libraries: eligible calls run on the GPU
// (Rblas_bridge_cuda.cpp, ../Rblas_cuda_common.cuh), the others on the CPU
// (../blas/Rsyrk_omp.cpp: lane-parallel OpenMP, else the reference code).
// All of them give the reference result bit for bit.

#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"
#include <cstdlib>

bool Rsyrk_omp(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *c, INTEGER const ldc);
void Rsyrk_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *c, INTEGER const ldc);
#if defined(MPLAPACK_BUILD_WITH_DD) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
// the 2010-2011 kernels (../../dd/cuda/Rsyrk_legacy_cuda.cu)
void Rsyrk_legacy(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *c, INTEGER const ldc);
#endif

void Rsyrk(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *c, INTEGER const ldc)
{
#if defined(MPLAPACK_BUILD_WITH_DD) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
    const char *legacy = std::getenv("MPLAPACK_DD_CUDA_LEGACY");
    if (legacy && legacy[0] == '1') {
        Rsyrk_legacy(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
        return;
    }
#endif
    if (MPLAPACK_CUDA_ENTRY(Rsyrk)(uplo, trans, n, k, alpha, a, lda, beta, c, ldc))
        return;
    if (Rsyrk_omp(uplo, trans, n, k, alpha, a, lda, beta, c, ldc))
        return;
    Rsyrk_reference_cpu(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}
