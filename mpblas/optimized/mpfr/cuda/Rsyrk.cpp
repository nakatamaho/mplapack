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

// Rsyrk of libmplapack_mpfr_opt_cuda: eligible calls run on the GPU
// (blas3_bridge_cuda.cpp), all others on the CPU
// (../../common/blas/Rsyrk_omp.cpp: lane-parallel OpenMP, else the reference
// code under the name Rsyrk_reference_cpu).

#include <mpblas_mpfr.h>

bool Rsyrk_mpfr_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class &beta, mpfr_class *c, mplapackint ldc);
void Rsyrk_reference_cpu(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);
bool Rsyrk_omp(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);

void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class const beta, mpfr_class *c, mplapackint const ldc)
{
    if (Rsyrk_mpfr_cuda(uplo, trans, n, k, alpha, a, lda, beta, c, ldc))
        return;
    if (Rsyrk_omp(uplo, trans, n, k, alpha, a, lda, beta, c, ldc))
        return;
    Rsyrk_reference_cpu(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}
