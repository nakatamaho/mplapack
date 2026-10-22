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

// The reference Ctbsv (mpblas/reference/Ctbsv.cpp, unchanged) under the name
// Ctbsv_reference_cpu, and Ctbsv_omp: for dd/td/qd, trans = N and incx = 1, the
// loops of the reference with the updates of x vectorized (they touch
// contiguous elements of the band column); the columns are solved one after
// the other, as in the reference, so the result is the reference one bit for
// bit.

#define Ctbsv Ctbsv_reference_cpu
#include "../../../reference/Ctbsv.cpp"
#undef Ctbsv

#include <algorithm>
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Cblas_simd_common.h"
#endif

bool Ctbsv_omp(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx)
{
#ifdef MPLAPACK_CBLAS_SIMD_AVAILABLE
    using namespace mplapack_blas;
    const bool upper = Mlsame(uplo, "U"), nounit = Mlsame(diag, "N");
    if (!Mlsame(trans, "N") || (!upper && !Mlsame(uplo, "L")) || (!nounit && !Mlsame(diag, "U")) || n <= 0 || k < 0 || lda < k + 1 || incx != 1)
        return false;
    if (blas_env_long("MPLAPACK_BLAS_SIMD", 1) == 0 || (double)n * (k + 1) < (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768))
        return false;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    for (long jj = 0; jj < n; jj++) {
        const long j = upper ? n - 1 - jj : jj;
        const COMPLEX *d = upper ? a + k + j * lda : a + j * lda;            // a(j, j)
        const long r0 = upper ? std::max(0L, j - (long)k) : j + 1;
        const long cnt = upper ? j - r0 : std::min((long)n - 1, j + (long)k) - j;
        const COMPLEX *col = upper ? a + (k + r0 - j) + j * lda : a + 1 + j * lda; // a(r0, j)
        if (x[j] != zero) {
            if (nounit)
                x[j] = x[j] / *d;
            const COMPLEX temp = x[j];
            cbk_axpy(cnt, x + r0, 1, temp, col, 1, true); // x(i) = x(i) - temp*a(i,j)
        }
    }
    return true;
#else
    (void)uplo, (void)trans, (void)diag, (void)n, (void)k, (void)a, (void)lda, (void)x, (void)incx;
    return false;
#endif
}
