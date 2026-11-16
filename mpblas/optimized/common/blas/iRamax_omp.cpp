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

#include "../Riamax_chunks_common.h"

// Goal 13: one pass over dx in fixed contiguous chunks (4096 elements), each
// chunk scanned serially, the chunk winners combined in chunk order.  The
// scores are the exact reference values and ties keep the first index, so the
// result is the reference index bit for bit for every thread count.  NaNs
// still go to the reference (unchanged contract).
bool iRamax_omp(INTEGER const n, REAL *dx, INTEGER const incx, INTEGER &result)
{
    using namespace mplapack_blas;
    if (n < 2 || incx != 1)
        return false;
    return iamax_chunks(n, [&](long i) { return abs(dx[i]); }, result);
}
