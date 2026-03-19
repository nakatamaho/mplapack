/*
 * Copyright (c) 2008-2025
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

// Derived from BLAS routine DNRM2.
// Original BLAS authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack_arithmetic_params.h>

REAL Rnrm2(INTEGER const n, REAL *x, INTEGER const incx) {
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
        ax = abs(x[ix - 1]);
        if (ax > tbig) {
            abig += pow2((ax * sbig));
            notbig = false;
        } else if (ax < tsml) {
            if (notbig) {
                asml += pow2((ax * ssml));
            }
        } else {
            amed += pow2(ax);
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
