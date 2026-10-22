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

// The reference Rnrm2 (mpblas/reference/Rnrm2.cpp, unchanged) under the name
// Rnrm2_reference_cpu, and Rnrm2_omp: the class (big, small, medium) and the term of every
// value in parallel, then the reference accumulation in the reference order
// (misc/gen_blas_l1_reductions.py).  Unit increment only.

#define Rnrm2 Rnrm2_reference_cpu
#include "../../../reference/Rnrm2.cpp"
#undef Rnrm2

#include <vector>
#include "../Rblas_split_common.h"

namespace {

// The reference Rnrm2 with its loop over x replaced by the accumulation of
// the precomputed classes and terms (the reference computes the term of a
// small value only while no big value has occurred; the term itself does
// not depend on that).
REAL Rnrm2_terms(INTEGER const n, REAL *x, INTEGER const incx, const signed char *cls, const REAL *term) {
    REAL return_value = 0.0;
    // Quick return if possible
    //
    const REAL zero = 0.0;
    return_value = zero;
    if (n <= 0) {
        return return_value;
    }
    //
    const REAL one = 1.0;
    REAL scl = one;
    REAL sumsq = zero;
    //
    // Compute the sum of squares in 3 accumulators:
    // abig -- sums of squares scaled down to avoid overflow
    // asml -- sums of squares scaled up to avoid underflow
    // amed -- sums of squares that do not require scaling
    // The thresholds and multipliers are
    // tbig -- values bigger than this are scaled down by sbig
    // tsml -- values smaller than this are scaled up by ssml
    //
    bool notbig = true;
    REAL asml = zero;
    REAL amed = zero;
    REAL abig = zero;
    INTEGER ix = 1;
    if (incx < 0) {
        ix = 1 - (n - 1) * incx;
    }
    INTEGER i = 0;
    REAL ax = 0.0;
    const auto &ap = mplapack::get_arithmetic_params<REAL>();
    const auto bp = mplapack::make_blue_scaling_params(ap);
    const REAL tbig = bp.tbig;
    const REAL sbig = bp.sbig;
    const REAL tsml = bp.tsml;
    const REAL ssml = bp.ssml;
    for (i = 1; i <= n; i = i + 1) {
        for (int c = 0; c < 1; c++) {
            const size_t t = (size_t)((ix - 1) * 1 + c);
            if (cls[t] == 0) {
                abig += term[t];
                notbig = false;
            } else if (cls[t] == 1) {
                if (notbig) {
                    asml += term[t];
                }
            } else {
                amed += term[t];
            }
        }
        ix += incx;
    }
    //
    // Combine abig and amed or amed and asml if more than one
    // accumulator was used.
    //
    const REAL maxn = ap.rmax;
    REAL ymin = 0.0;
    REAL ymax = 0.0;
    if (abig > zero) {
        //
        // Combine abig and amed if abig > 0.
        //
        if ((amed > zero) || (amed > maxn) || (amed != amed)) {
            abig += (amed * sbig) * sbig;
        }
        scl = one / sbig;
        sumsq = abig;
    } else if (asml > zero) {
        //
        // Combine amed and asml if asml > 0.
        //
        if ((amed > zero) || (amed > maxn) || (amed != amed)) {
            amed = sqrt(amed);
            asml = sqrt(asml) / ssml;
            if (asml > amed) {
                ymin = amed;
                ymax = asml;
            } else {
                ymin = asml;
                ymax = amed;
            }
            scl = one;
            sumsq = pow2(ymax) * (one + pow2((ymin / ymax)));
        } else {
            scl = one / ssml;
            sumsq = asml;
        }
    } else {
        //
        // Otherwise all values are mid-range
        //
        scl = one;
        sumsq = amed;
    }
    return_value = scl * sqrt(sumsq);
    return return_value;
}
} // namespace

bool Rnrm2_omp(INTEGER const n, REAL *x, INTEGER const incx, REAL &result)
{
    using namespace mplapack_blas;
    if (n <= 0 || incx != 1)
        return false;
    std::vector<signed char> cls((size_t)n * 1);
    std::vector<REAL> term((size_t)n * 1);
    const auto &ap = mplapack::get_arithmetic_params<REAL>();
    const auto bp = mplapack::make_blue_scaling_params(ap);
    const REAL tbig = bp.tbig;
    const REAL sbig = bp.sbig;
    const REAL tsml = bp.tsml;
    const REAL ssml = bp.ssml;
    if (!blas_split(4.0 * n * 1, n, [&](long k0, long k1) {
            REAL ax = 0.0;
            for (long k = k0; k < k1; k++) {
                size_t t = (size_t)k * 1;
                ax = abs(x[k]);
                if (ax > tbig) {
                    cls[t] = 0, term[t] = pow2((ax * sbig));
                } else if (ax < tsml) {
                    cls[t] = 1, term[t] = pow2((ax * ssml));
                } else {
                    cls[t] = 2, term[t] = pow2(ax);
                }
                t++;
            }
        }))
        return false;
    result = Rnrm2_terms(n, x, incx, cls.data(), term.data());
    return true;
}
