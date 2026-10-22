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

// The reference Rtpmv (mpblas/reference/Rtpmv.cpp, unchanged) under the name
// Rtpmv_reference_cpu, and Rtpmv_omp: ranges of x in parallel, from a copy of the original x.
// Every element receives the operations of the reference code in the same
// order, so the result is the reference one bit for bit.  Unit increments
// only; the other calls run the reference code.

#define Rtpmv Rtpmv_reference_cpu
#include "../../../reference/Rtpmv.cpp"
#undef Rtpmv

#include <algorithm>
#include <vector>
#include "../Rblas_split_common.h"

namespace {

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline INTEGER imax(INTEGER a, INTEGER b) { return a > b ? a : b; }
inline INTEGER imin(INTEGER a, INTEGER b) { return a < b ? a : b; }

// Rtpmv restricted to the elements r0 < i <= r1 of x (1-based); x0 is the
// original x (the reference reads x(j) before any update of it, so it is x0(j)).
void Rtpmv_rows(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *ap, REAL *x, INTEGER const incx, const REAL *x0, INTEGER const r0, INTEGER const r1)
{
    bool nounit = Mlsame(diag, "N");
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL temp = 0.0;
    INTEGER i = 0;
    if (Mlsame(trans, "N")) {
        if (Mlsame(uplo, "U")) {
            for (j = r0 + 1; j <= n; j = j + 1) {
                if (x0[j - 1] != zero) {
                    temp = x0[j - 1];
                    for (i = imax(1, r0 + 1); i <= imin(j - 1, r1); i = i + 1) {
                        x[i - 1] += temp * ap[(j * (j - 1) / 2 + i) - 1];
                    }
                    if (nounit && j <= r1) {
                        x[j - 1] = x[j - 1] * ap[(j * (j - 1) / 2 + j) - 1];
                    }
                }
            }
        } else {
            for (j = r1; j >= 1; j = j - 1) {
                if (x0[j - 1] != zero) {
                    temp = x0[j - 1];
                    for (i = imin(n, r1); i >= imax(j + 1, r0 + 1); i = i - 1) {
                        x[i - 1] += temp * ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + i - j + 1) - 1];
                    }
                    if (nounit && j > r0) {
                        x[j - 1] = x[j - 1] * ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + j - j + 1) - 1];
                    }
                }
            }
        }
    } else {
        // x(j) from column j and the original x
        if (Mlsame(uplo, "U")) {
            for (j = r1; j >= r0 + 1; j = j - 1) {
                temp = x0[j - 1];
                if (nounit) {
                    temp = temp * ap[(j * (j - 1) / 2 + j) - 1];
                }
                for (i = j - 1; i >= 1; i = i - 1) {
                    temp += ap[(j * (j - 1) / 2 + i) - 1] * x0[i - 1];
                }
                x[j - 1] = temp;
            }
        } else {
            for (j = r0 + 1; j <= r1; j = j + 1) {
                temp = x0[j - 1];
                if (nounit) {
                    temp = temp * ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + j - j + 1) - 1];
                }
                for (i = j + 1; i <= n; i = i + 1) {
                    temp += ap[((j - 1) * n - (j - 1) * (j - 2) / 2 + i - j + 1) - 1] * x0[i - 1];
                }
                x[j - 1] = temp;
            }
        }
    }
}

} // namespace

bool Rtpmv_omp(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *ap, REAL *x, INTEGER const incx)
{
    using namespace mplapack_blas;
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1)
        return false;
    const std::vector<REAL> x0(x, x + n);
    return blas_split(4.0 * n * n, n, [&](long p0, long p1) { Rtpmv_rows(uplo, trans, diag, n, ap, x, incx, x0.data(), (INTEGER)p0, (INTEGER)p1); });
}
