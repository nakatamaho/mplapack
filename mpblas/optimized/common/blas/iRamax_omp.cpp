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

// The reference iRamax (mpblas/reference/iRamax.cpp, unchanged) under the name
// iRamax_reference_cpu, and iRamax_omp: the values in parallel; the first index of their
// maximum, as the reference finds it (comparisons are exact; values that
// compare false with everything, NaN, go to the reference)
// (misc/gen_blas_l1_reductions.py).  Unit increment only.

#define iRamax iRamax_reference_cpu
#include "../../../reference/iRamax.cpp"
#undef iRamax

#include <vector>
#include "../Rblas_split_common.h"

bool iRamax_omp(INTEGER const n, REAL *dx, INTEGER const incx, INTEGER &result)
{
    using namespace mplapack_blas;
    if (n < 2 || incx != 1)
        return false;
    std::vector<REAL> v((size_t)n);
    std::vector<char> ok((size_t)n);
    if (!blas_split(2.0 * n, n, [&](long i0, long i1) {
            for (long i = i0; i < i1; i++) {
                v[(size_t)i] = abs(dx[i]);
                ok[(size_t)i] = v[(size_t)i] == v[(size_t)i];
            }
        }))
        return false;
    for (INTEGER i = 0; i < n; i++)
        if (!ok[(size_t)i])
            return false;
    // the reference: dmax = v(1); for i = 2..n: if (v(i) > dmax) take i
    INTEGER best = 1;
    REAL dmax = v[0];
    for (INTEGER i = 2; i <= n; i = i + 1) {
        if (v[(size_t)(i - 1)] > dmax) {
            best = i;
            dmax = v[(size_t)(i - 1)];
        }
    }
    result = best;
    return true;
}
