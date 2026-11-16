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

// The reference Cgeru (mpblas/reference/Cgeru.cpp, unchanged) under the name
// Cgeru_reference_cpu, and Cgeru_omp: the reference code on independent parts of
// the output in parallel.

#define Cgeru Cgeru_reference_cpu
#include "../../../reference/Cgeru.cpp"
#undef Cgeru

#include <algorithm>
#include "../Rblas_split_common.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Cblas_simd_common.h"
#endif

bool Cgeru_omp(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    using namespace mplapack_blas;
    if (m <= 0 || n <= 0 || incx == 0 || incy <= 0 || lda < std::max((INTEGER)1, m))
        return false;
#ifdef MPLAPACK_CBLAS_SIMD_AVAILABLE
    if (cger_simd<COMPLEX>(false, m, n, alpha, x, incx, y, incy, a, lda))
        return true;
#endif
    return blas_split(4.0 * m * n, n, [&](long p0, long p1) { Cgeru_reference_cpu(m, p1 - p0, alpha, x, incx, y + p0 * incy, incy, a + p0 * lda, lda); });
}
