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

// Rlaswp of the optimized libraries: the row interchanges of the columns are
// independent, so large calls run column ranges in parallel, each with the
// reference code (mplapack/reference/Rlaswp.cpp, unchanged, under the name
// Rlaswp_reference_cpu); the result is the reference result bit for bit.

#define Rlaswp Rlaswp_reference_cpu
#include "../../reference/Rlaswp.cpp"
#undef Rlaswp

#include "../../../mpblas/optimized/common/Rblas_split_common.h"

void Rlaswp(INTEGER const n, REAL *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx)
{
    if (n > 0 && incx != 0 && k2 >= k1 && mplapack_blas::blas_split(4.0 * n * (k2 - k1 + 1), n, [&](long p0, long p1) { Rlaswp_reference_cpu(p1 - p0, a + p0 * lda, lda, k1, k2, ipiv, incx); }))
        return;
    Rlaswp_reference_cpu(n, a, lda, k1, k2, ipiv, incx);
}
