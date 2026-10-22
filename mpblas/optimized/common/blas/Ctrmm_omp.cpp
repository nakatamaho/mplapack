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

// The reference Ctrmm (mpblas/reference/Ctrmm.cpp, unchanged) under the name
// Ctrmm_reference_cpu, and Ctrmm_omp: the reference code on independent parts of
// the output in parallel.

#define Ctrmm Ctrmm_reference_cpu
#include "../../../reference/Ctrmm.cpp"
#undef Ctrmm

#include <algorithm>
#include "../Rblas_split_common.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Cblas_simd_common.h"
#endif

namespace {
bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
} // namespace

bool Ctrmm_omp(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb)
{
    using namespace mplapack_blas;
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(transa, "N") && !ls(transa, "T") && !ls(transa, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || m <= 0 || n <= 0)
        return false;
    if (lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m))
        return false;
#ifdef MPLAPACK_CBLAS_SIMD_AVAILABLE
    if (ctrxm_simd<COMPLEX>(false, left, ls(uplo, "U"), ls(transa, "N"), ls(transa, "T"), ls(diag, "N"), m, n, alpha, a, lda, b, ldb))
        return true;
#endif
    const double work = 4.0 * m * n * (left ? m : n);
    if (left) // columns of B
        return blas_split(work, n, [&](long p0, long p1) { Ctrmm_reference_cpu(side, uplo, transa, diag, m, p1 - p0, alpha, a, lda, b + p0 * ldb, ldb); });
    return blas_split(work, m, [&](long p0, long p1) { Ctrmm_reference_cpu(side, uplo, transa, diag, p1 - p0, n, alpha, a, lda, b + p0, ldb); });
}
