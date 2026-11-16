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

// The reference Rsyrk (mpblas/reference/Rsyrk.cpp, unchanged) under the name
// Rsyrk_reference_cpu, and Rsyrk_omp: the same operations in parallel lanes.

#define Rsyrk Rsyrk_reference_cpu
#include "../../../reference/Rsyrk.cpp"
#undef Rsyrk

#include "../Rblas_cpu_common.h"

namespace {
bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
} // namespace

bool Rsyrk_omp(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *cc, INTEGER const ldc)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_syrk<REAL>(c, ls, false, uplo, trans, n, k, alpha, a, lda, (REAL *)0, 1, beta, cc, ldc))
        return false;
    return blas_cpu(c);
}
