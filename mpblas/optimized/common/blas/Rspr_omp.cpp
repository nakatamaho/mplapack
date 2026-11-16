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

// The reference Rspr (mpblas/reference/Rspr.cpp, unchanged) under the name
// Rspr_reference_cpu, and Rspr_omp: ranges of columns of AP in parallel.
// Every element receives the operations of the reference code in the same
// order, so the result is the reference one bit for bit.  Unit increments
// only; the other calls run the reference code.

#define Rspr Rspr_reference_cpu
#include "../../../reference/Rspr.cpp"
#undef Rspr

#include <algorithm>
#include <vector>
#include "../Rblas_split_common.h"

namespace {

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline INTEGER imax(INTEGER a, INTEGER b) { return a > b ? a : b; }
inline INTEGER imin(INTEGER a, INTEGER b) { return a < b ? a : b; }

// Rspr restricted to the columns r0 < j <= r1 of AP (1-based); kk is the
// position of the first element of column j, as the reference counts it.
void Rspr_rows(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *ap, INTEGER const r0, INTEGER const r1)
{
    const REAL zero = 0.0;
    INTEGER kk = 1;
    INTEGER j = 0;
    REAL temp = 0.0;
    INTEGER k = 0;
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        for (j = r0 + 1; j <= r1; j = j + 1) {
            kk = j * (j - 1) / 2 + 1;
            if (x[j - 1] != zero) {
                temp = alpha * x[j - 1];
                k = kk;
                for (i = 1; i <= j; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp;
                    k++;
                }
            }
        }
    } else {
        for (j = r0 + 1; j <= r1; j = j + 1) {
            kk = (j - 1) * n - (j - 1) * (j - 2) / 2 + 1;
            if (x[j - 1] != zero) {
                temp = alpha * x[j - 1];
                k = kk;
                for (i = j; i <= n; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp;
                    k++;
                }
            }
        }
    }
}

} // namespace

bool Rspr_omp(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *ap)
{
    using namespace mplapack_blas;
    const REAL zero = 0.0;
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || alpha == zero || incx != 1)
        return false;
    return blas_split(2.0 * n * n, n, [&](long p0, long p1) { Rspr_rows(uplo, n, alpha, x, incx, ap, (INTEGER)p0, (INTEGER)p1); });
}
