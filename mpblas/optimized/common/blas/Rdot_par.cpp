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

// The reference Rdot (mpblas/reference/Rdot.cpp, unchanged) under the name
// Rdot_reference_cpu, and Rdot_par: for unit increments, the sums of five of
// the unrolled reference loop in parallel, then their sum in the reference
// order, so the result is the reference one bit for bit.  Called by the
// Rdot dispatchers of the optimized libraries (<type>/Rdot.cpp).

#define Rdot Rdot_reference_cpu
#include "../../../reference/Rdot.cpp"
#undef Rdot

#include <vector>
#include "../Rblas_split_common.h"

bool Rdot_par(INTEGER const n, REAL *dx, INTEGER const incx, REAL *dy, INTEGER const incy, REAL &result)
{
    using namespace mplapack_blas;
    if (n < 10 || incx != 1 || incy != 1)
        return false;
    const INTEGER m = n % 5, groups = (n - m) / 5;
    std::vector<REAL> g((size_t)groups);
    if (!blas_split(4.0 * n, groups, [&](long q0, long q1) {
            for (long q = q0; q < q1; q++) {
                const INTEGER i = m + 1 + 5 * q;
                g[(size_t)q] = dx[i - 1] * dy[i - 1] + dx[(i + 1) - 1] * dy[(i + 1) - 1] + dx[(i + 2) - 1] * dy[(i + 2) - 1] + dx[(i + 3) - 1] * dy[(i + 3) - 1] + dx[(i + 4) - 1] * dy[(i + 4) - 1];
            }
        }))
        return false;
    REAL dtemp = 0.0;
    for (INTEGER i = 1; i <= m; i = i + 1) {
        dtemp += dx[i - 1] * dy[i - 1];
    }
    for (INTEGER q = 0; q < groups; q++) {
        dtemp += g[(size_t)q];
    }
    result = dtemp;
    return true;
}
