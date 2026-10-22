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

// The reference Ctrsv (mpblas/reference/Ctrsv.cpp, unchanged) under the name
// Ctrsv_reference_cpu, and Ctrsv_omp: for dd/td/qd and incx = 1 the loops of the
// reference with AVX2 (misc/gen_blas_tri_solves.py): trans = N, the updates
// of x in a column; trans = T, C, the products of a column, then
// subtracted in the reference order.  The columns are solved one after the
// other, as in the reference.

#define Ctrsv Ctrsv_reference_cpu
#include "../../../reference/Ctrsv.cpp"
#undef Ctrsv

#include <algorithm>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Cblas_simd_common.h"
#endif

bool Ctrsv_omp(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx)
{
#ifdef MPLAPACK_CBLAS_SIMD_AVAILABLE
    using namespace mplapack_blas;
    const bool upper = Mlsame(uplo, "U"), nounit = Mlsame(diag, "N"), notr = Mlsame(trans, "N");
    if (!notr && !Mlsame(trans, "T") && !Mlsame(trans, "C"))
        return false;
    if ((!upper && !Mlsame(uplo, "L")) || (!nounit && !Mlsame(diag, "U")) || n <= 0 || lda < std::max((INTEGER)1, n) || incx != 1)
        return false;
    if (blas_env_long("MPLAPACK_BLAS_SIMD", 1) == 0 || (double)n * n < (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768))
        return false;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    const bool noconj = Mlsame(trans, "T");
    if (notr) {
        for (long jj = 0; jj < n; jj++) {
            const long j = upper ? n - 1 - jj : jj;
            const COMPLEX *d = a + j + j * lda; // a(j, j)
            const long r0 = upper ? 0 : j + 1;
            const long cnt = upper ? j : (long)n - j - 1;
            const COMPLEX *col = a + r0 + j * lda; // a(r0, j)
            if (x[j] != zero) {
                if (nounit)
                    x[j] = x[j] / *d;
                const COMPLEX temp = x[j];
                cbk_axpy(cnt, x + r0, 1, temp, col, 1, true); // x(i) = x(i) - temp*a(i,j)
            }
        }
        return true;
    }
    // x(j) = (x(j) - a(r0, j)*x(r0) - ...) / a(j, j): upper i ascending, lower
    // i descending, as the reference
    std::vector<COMPLEX> p((size_t)n);
    for (long jj = 0; jj < n; jj++) {
        const long j = upper ? jj : n - 1 - jj;
        const COMPLEX *d = a + j + j * lda; // a(j, j)
        const long r0 = upper ? 0 : j + 1;
        const long cnt = upper ? j : (long)n - j - 1;
        const COMPLEX *col = a + r0 + j * lda; // a(r0, j)
        cbk_vmul(cnt, p.data(), col, x + r0, !noconj);
        COMPLEX temp = x[j];
        if (upper)
            for (long t = 0; t < cnt; t++)
                temp = temp - p[(size_t)t];
        else
            for (long t = cnt - 1; t >= 0; t--)
                temp = temp - p[(size_t)t];
        if (nounit)
            temp = noconj ? temp / *d : temp / conj(*d);
        x[j] = temp;
    }
    return true;
#else
    (void)uplo, (void)trans, (void)diag, (void)n, (void)a, (void)lda, (void)x, (void)incx;
    return false;
#endif
}
