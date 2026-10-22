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

// The reference Rtrsv (mpblas/reference/Rtrsv.cpp, unchanged) under the name
// Rtrsv_reference_cpu, and Rtrsv_omp: for dd/td/qd, trans = N and incx = 1, the
// loops of the reference with the updates of x vectorized; the columns are
// solved one after the other, as in the reference.

#define Rtrsv Rtrsv_reference_cpu
#include "../../../reference/Rtrsv.cpp"
#undef Rtrsv

#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Rblas_simd_common.h"
#endif

bool Rtrsv_omp(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx)
{
#ifdef MPLAPACK_BLAS_SIMD_AVAILABLE
    using namespace mplapack_blas;
    const bool upper = Mlsame(uplo, "U"), nounit = Mlsame(diag, "N");
    if (!Mlsame(trans, "N") || (!upper && !Mlsame(uplo, "L")) || (!nounit && !Mlsame(diag, "U")) || n <= 0 || lda < std::max((INTEGER)1, n) || incx != 1)
        return false;
    if (blas_env_long("MPLAPACK_BLAS_SIMD", 1) == 0 || (double)n * n < (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768))
        return false;
    const REAL zero = 0.0;
    for (long jj = 0; jj < n; jj++) {
        const long j = upper ? n - 1 - jj : jj;
        if (x[j] != zero) {
            if (nounit)
                x[j] = x[j] / a[j + j * lda];
            const REAL temp = x[j];
            const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1;
            bk_axpy(cnt, x + r0, 1, temp, a + r0 + j * lda, 1, true); // x(i) = x(i) - temp*a(i,j)
        }
    }
    return true;
#else
    (void)uplo, (void)trans, (void)diag, (void)n, (void)a, (void)lda, (void)x, (void)incx;
    return false;
#endif
}
