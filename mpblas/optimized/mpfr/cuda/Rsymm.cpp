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

// Rsymm of libmplapack_mpfr_opt_cuda: eligible calls run on the GPU
// (blas3_bridge_cuda.cpp), all others in the reference implementation
// (Rsymm_ref_cuda.cpp: ../../../reference/Rsymm.cpp under another name).

#include <mpblas_mpfr.h>

bool Rsymm_mpfr_cuda(const char *side, const char *uplo, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class *b, mplapackint ldb, const mpfr_class &beta, mpfr_class *c, mplapackint ldc);
void Rsymm_reference_cpu(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);

void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc)
{
    if (Rsymm_mpfr_cuda(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc))
        return;
    Rsymm_reference_cpu(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}
