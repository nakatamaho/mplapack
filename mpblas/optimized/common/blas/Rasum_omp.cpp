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

// The reference Rasum (mpblas/reference/Rasum.cpp, unchanged) under the name
// Rasum_reference_cpu, and Rasum_omp: the sums of six of the unrolled reference loop in
// parallel, then their sum in the reference order
// (misc/gen_blas_l1_reductions.py).  Unit increment only.

#define Rasum Rasum_reference_cpu
#include "../../../reference/Rasum.cpp"
#undef Rasum

#include <vector>
#include "../Rblas_split_common.h"

bool Rasum_omp(INTEGER const n, REAL *dx, INTEGER const incx, REAL &result)
{
    using namespace mplapack_blas;
    if (n < 12 || incx != 1)
        return false;
    const INTEGER m = n % 6, groups = (n - m) / 6;
    std::vector<REAL> g((size_t)groups);
    if (!blas_split(2.0 * n, groups, [&](long q0, long q1) {
            for (long q = q0; q < q1; q++) {
                const INTEGER i = m + 1 + 6 * q;
                g[(size_t)q] = abs(dx[i - 1]) + abs(dx[(i + 1) - 1]) + abs(dx[(i + 2) - 1]) + abs(dx[(i + 3) - 1]) + abs(dx[(i + 4) - 1]) + abs(dx[(i + 5) - 1]);
            }
        }))
        return false;
    REAL dtemp = 0.0;
    for (INTEGER i = 1; i <= m; i = i + 1) {
        dtemp += abs(dx[i - 1]);
    }
    for (INTEGER q = 0; q < groups; q++) {
        dtemp += g[(size_t)q];
    }
    result = dtemp;
    return true;
}
