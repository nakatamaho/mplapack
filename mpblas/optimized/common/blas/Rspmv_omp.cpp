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

// The reference Rspmv (mpblas/reference/Rspmv.cpp, unchanged) under the name
// Rspmv_reference_cpu, and Rspmv_omp: ranges of rows of y in parallel.
// Every element receives the operations of the reference code in the same
// order, so the result is the reference one bit for bit.  Unit increments
// only; the other calls run the reference code.

#define Rspmv Rspmv_reference_cpu
#include "../../../reference/Rspmv.cpp"
#undef Rspmv

#include <algorithm>
#include <vector>
#include "../Rblas_split_common.h"

namespace {

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline INTEGER imax(INTEGER a, INTEGER b) { return a > b ? a : b; }
inline INTEGER imin(INTEGER a, INTEGER b) { return a < b ? a : b; }

// Rspmv restricted to the elements r0 < i <= r1 of y (1-based).
void Rspmv_rows(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy, INTEGER const r0, INTEGER const r1)
{
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER i = 0;
    if (beta != one) {
        if (beta == zero) {
            for (i = r0 + 1; i <= r1; i = i + 1) {
                y[i - 1] = zero;
            }
        } else {
            for (i = r0 + 1; i <= r1; i = i + 1) {
                y[i - 1] = beta * y[i - 1];
            }
        }
    }
    INTEGER j = 0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    if (Mlsame(uplo, "U")) {
        // row i: column i (with the whole column above i), then columns j > i
        for (j = r0 + 1; j <= n; j = j + 1) {
            temp1 = alpha * x[j - 1];
            for (i = imax(1, r0 + 1); i <= imin(j - 1, r1); i = i + 1) {
                y[i - 1] += temp1 * ap[(j * (j - 1) / 2 + i) - 1];
            }
            if (j <= r1) {
                temp2 = zero;
                for (i = 1; i <= j - 1; i = i + 1) {
                    temp2 += ap[(j * (j - 1) / 2 + i) - 1] * x[i - 1];
                }
                y[j - 1] += temp1 * ap[(j * (j - 1) / 2 + j) - 1] + alpha * temp2;
            }
        }
    } else {
        // row i: columns j < i, then column i (with the whole column below i)
        for (j = 1; j <= r1; j = j + 1) {
            temp1 = alpha * x[j - 1];
            if (j > r0) {
                y[j - 1] += temp1 * ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + j - j + 1) - 1];
            }
            for (i = imax(j + 1, r0 + 1); i <= imin(n, r1); i = i + 1) {
                y[i - 1] += temp1 * ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + i - j + 1) - 1];
            }
            if (j > r0) {
                temp2 = zero;
                for (i = j + 1; i <= n; i = i + 1) {
                    temp2 += ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + i - j + 1) - 1] * x[i - 1];
                }
                y[j - 1] += alpha * temp2;
            }
        }
    }
}

} // namespace

bool Rspmv_omp(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    using namespace mplapack_blas;
    const REAL zero = 0.0;
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == zero) // alpha == 0: only y = beta y
        return false;
    return blas_split(4.0 * n * n, n, [&](long p0, long p1) { Rspmv_rows(uplo, n, alpha, ap, x, incx, beta, y, incy, (INTEGER)p0, (INTEGER)p1); });
}
