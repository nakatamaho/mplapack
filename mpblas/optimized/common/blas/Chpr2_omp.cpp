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

// The reference Chpr2 (mpblas/reference/Chpr2.cpp, unchanged) under the name
// Chpr2_reference_cpu, and Chpr2_omp: ranges of columns of AP in parallel.
// Every element receives the operations of the reference code in the same
// order, so the result is the reference one bit for bit.  Unit increments
// only; the other calls run the reference code.

#define Chpr2 Chpr2_reference_cpu
#include "../../../reference/Chpr2.cpp"
#undef Chpr2

#include <algorithm>
#include <vector>
#include "../Rblas_split_common.h"

namespace {

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline INTEGER imax(INTEGER a, INTEGER b) { return a > b ? a : b; }
inline INTEGER imin(INTEGER a, INTEGER b) { return a < b ? a : b; }

// Chpr2 restricted to the columns r0 < j <= r1 of AP (1-based); kk is the
// position of the first element of column j, as the reference counts it.
void Chpr2_rows(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *ap, INTEGER const r0, INTEGER const r1)
{
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    INTEGER kk = 1;
    INTEGER j = 0;
    COMPLEX temp1 = 0.0;
    COMPLEX temp2 = 0.0;
    INTEGER k = 0;
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        for (j = r0 + 1; j <= r1; j = j + 1) {
            kk = j * (j - 1) / 2 + 1;
            if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                temp1 = alpha * conj(y[j - 1]);
                temp2 = conj(alpha * x[j - 1]);
                k = kk;
                for (i = 1; i <= j - 1; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                    k++;
                }
                ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1].real() + (x[j - 1] * temp1 + y[j - 1] * temp2).real();
            } else {
                ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1].real();
            }
        }
    } else {
        for (j = r0 + 1; j <= r1; j = j + 1) {
            kk = (j - 1) * n - (j - 1) * (j - 2) / 2 + 1;
            if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                temp1 = alpha * conj(y[j - 1]);
                temp2 = conj(alpha * x[j - 1]);
                ap[kk - 1] = ap[kk - 1].real() + (x[j - 1] * temp1 + y[j - 1] * temp2).real();
                k = kk + 1;
                for (i = j + 1; i <= n; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                    k++;
                }
            } else {
                ap[kk - 1] = ap[kk - 1].real();
            }
        }
    }
}

} // namespace

bool Chpr2_omp(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *ap)
{
    using namespace mplapack_blas;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || alpha == zero || incx != 1 || incy != 1)
        return false;
    return blas_split(4.0 * n * n, n, [&](long p0, long p1) { Chpr2_rows(uplo, n, alpha, x, incx, y, incy, ap, (INTEGER)p0, (INTEGER)p1); });
}
