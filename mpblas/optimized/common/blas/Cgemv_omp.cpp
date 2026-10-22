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

// The reference Cgemv (mpblas/reference/Cgemv.cpp, unchanged) under the name
// Cgemv_reference_cpu, and Cgemv_omp: the reference code on independent parts of
// the output in parallel.

#define Cgemv Cgemv_reference_cpu
#include "../../../reference/Cgemv.cpp"
#undef Cgemv

#include <algorithm>
#include "../Rblas_split_common.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Cblas_simd_common.h"
#endif

namespace {
bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
} // namespace

bool Cgemv_omp(const char *trans, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    using namespace mplapack_blas;
    const bool notr = ls(trans, "N");
    if ((!notr && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, m) || incx == 0 || incy <= 0)
        return false;
#ifdef MPLAPACK_CBLAS_SIMD_AVAILABLE
    if (notr && cgemv_simd<COMPLEX>(m, n, alpha, a, lda, x, incx, beta, y, incy))
        return true;
#endif
    const double work = 4.0 * m * n;
    if (notr) // y(i), i < m
        return blas_split(work, m, [&](long p0, long p1) { Cgemv_reference_cpu(trans, p1 - p0, n, alpha, a + p0, lda, x, incx, beta, y + p0 * incy, incy); });
    return blas_split(work, n, [&](long p0, long p1) { Cgemv_reference_cpu(trans, m, p1 - p0, alpha, a + p0 * lda, lda, x, incx, beta, y + p0 * incy, incy); });
}
