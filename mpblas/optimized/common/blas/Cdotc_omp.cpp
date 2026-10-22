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

// The reference Cdotc (mpblas/reference/Cdotc.cpp, unchanged) under the name
// Cdotc_reference_cpu, and Cdotc_omp: the products in parallel, then the sum in the
// reference order
// (misc/gen_blas_l1_reductions.py).  Unit increment only.

#define Cdotc Cdotc_reference_cpu
#include "../../../reference/Cdotc.cpp"
#undef Cdotc

#include <vector>
#include "../Rblas_split_common.h"

bool Cdotc_omp(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy, COMPLEX &result)
{
    using namespace mplapack_blas;
    if (n <= 0 || incx != 1 || incy != 1)
        return false;
    std::vector<COMPLEX> p((size_t)n);
    if (!blas_split(8.0 * n, n, [&](long i0, long i1) {
            for (long i = i0; i < i1; i++)
                p[(size_t)i] = conj(zx[i]) * zy[i];
        }))
        return false;
    COMPLEX ztemp = COMPLEX(0.0, 0.0);
    for (INTEGER i = 1; i <= n; i = i + 1) {
        ztemp += p[(size_t)(i - 1)];
    }
    result = ztemp;
    return true;
}
