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

// The reference Chemm (mpblas/reference/Chemm.cpp, unchanged) under the name
// Chemm_reference_cpu, and Chemm_omp: the reference code on independent parts of
// the output in parallel.

#define Chemm Chemm_reference_cpu
#include "../../../reference/Chemm.cpp"
#undef Chemm

#include <algorithm>
#include "../Rblas_split_common.h"
#include "../Rblas_cpu_common.h"

namespace {
bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
} // namespace

bool Chemm_omp(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    using namespace mplapack_blas;
    const bool left = ls(side, "L"), upper = ls(uplo, "U");
    if ((!left && !ls(side, "R")) || (!upper && !ls(uplo, "L")) || m <= 0 || n <= 0)
        return false;
    if (lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m) || ldc < std::max((INTEGER)1, m))
        return false;
    // the small-block engine where it is faster (../Rblas_block_common.h)
    if (!(alpha == COMPLEX(0.0, 0.0)) && blk_symm<blas_thread_context, COMPLEX, true>(left, upper, m, n, alpha, a, lda, b, ldb, beta, c, ldc))
        return true;
    const double work = 4.0 * m * n * (left ? m : n);
    if (left) // columns of C
        return blas_split(work, n, [&](long p0, long p1) { Chemm_reference_cpu(side, uplo, m, p1 - p0, alpha, a, lda, b + p0 * ldb, ldb, beta, c + p0 * ldc, ldc); });
    return blas_split(work, m, [&](long p0, long p1) { Chemm_reference_cpu(side, uplo, p1 - p0, n, alpha, a, lda, b + p0, ldb, beta, c + p0, ldc); });
}
