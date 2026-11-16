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

// The reference Cgbmv (mpblas/reference/Cgbmv.cpp, unchanged) under the name
// Cgbmv_reference_cpu, and Cgbmv_omp: ranges of rows of y in parallel.
// Every element receives the operations of the reference code in the same
// order, so the result is the reference one bit for bit.  Unit increments
// only; the other calls run the reference code.

#define Cgbmv Cgbmv_reference_cpu
#include "../../../reference/Cgbmv.cpp"
#undef Cgbmv

#include <algorithm>
#include <vector>
#include "../Rblas_split_common.h"

namespace {

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline INTEGER imax(INTEGER a, INTEGER b) { return a > b ? a : b; }
inline INTEGER imin(INTEGER a, INTEGER b) { return a < b ? a : b; }

// Cgbmv restricted to the elements r0 < i <= r1 of y (1-based).
void Cgbmv_rows(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy, INTEGER const r0, INTEGER const r1)
{
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
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
    bool noconj = Mlsame(trans, "T");
    const INTEGER kup1 = ku + 1;
    INTEGER j = 0;
    COMPLEX temp = 0.0;
    INTEGER k = 0;
    if (Mlsame(trans, "N")) {
        // row i: the columns j of its band, in order
        for (j = imax((INTEGER)1, r0 + 1 - kl); j <= imin(n, r1 + ku); j = j + 1) {
            temp = alpha * x[j - 1];
            k = kup1 - j;
            for (i = imax(imax((INTEGER)1, j - ku), r0 + 1); i <= imin(imin(m, j + kl), r1); i = i + 1) {
                y[i - 1] += temp * a[((k + i) - 1) + (j - 1) * lda];
            }
        }
    } else {
        // y(j): column j
        for (j = r0 + 1; j <= r1; j = j + 1) {
            temp = zero;
            k = kup1 - j;
            if (noconj) {
                for (i = imax((INTEGER)1, j - ku); i <= imin(m, j + kl); i = i + 1) {
                    temp += a[((k + i) - 1) + (j - 1) * lda] * x[i - 1];
                }
            } else {
                for (i = imax((INTEGER)1, j - ku); i <= imin(m, j + kl); i = i + 1) {
                    temp += conj(a[((k + i) - 1) + (j - 1) * lda]) * x[i - 1];
                }
            }
            y[j - 1] += alpha * temp;
        }
    }
}

} // namespace

bool Cgbmv_omp(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    using namespace mplapack_blas;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    if ((!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || kl < 0 || ku < 0 || lda < kl + ku + 1 || incx != 1 || incy != 1 || alpha == zero) // alpha == 0: only y = beta y
        return false;
    return blas_split(4.0 * (ls(trans, "N") ? m : n) * (kl + ku + 1), ls(trans, "N") ? m : n, [&](long p0, long p1) { Cgbmv_rows(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, (INTEGER)p0, (INTEGER)p1); });
}
