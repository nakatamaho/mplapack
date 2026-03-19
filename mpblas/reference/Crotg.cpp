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

// Derived from BLAS routine ZROTG.
// Original BLAS authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack_arithmetic_params.h>

inline REAL abssq(COMPLEX ff) {
    REAL temp;
    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
    return temp;
}

void Crotg(COMPLEX &a, COMPLEX const b, REAL &c, COMPLEX &s) {
    //
    COMPLEX t = 0.0;
    //
    COMPLEX f = a;
    COMPLEX g = b;
    const COMPLEX czero = 0.0;
    const REAL one = 1.0;
    COMPLEX r = 0.0;
    const REAL zero = 0.0;
    REAL g1 = 0.0;
    const auto &ap = mplapack::get_arithmetic_params<REAL>();
    const REAL safmin = ap.safmin;
    REAL rtmax = 0.0;
    const REAL safmax = ap.safmax;
    const REAL rtmin = sqrt(safmin);
    REAL g2 = 0.0;
    REAL d = 0.0;
    REAL u = 0.0;
    COMPLEX gs = 0.0;
    REAL f1 = 0.0;
    REAL f2 = 0.0;
    REAL h2 = 0.0;
    REAL v = 0.0;
    REAL w = 0.0;
    COMPLEX fs = 0.0;
    if (g == czero) {
        c = one;
        s = czero;
        r = f;
    } else if (f == czero) {
        c = zero;
        if (g.real() == zero) {
            r = abs(g.imag());
            s = conj(g) / r;
        } else if (g.imag() == zero) {
            r = abs(g.real());
            s = conj(g) / r;
        } else {
            g1 = max(abs(g.real()), abs(g.imag()));
            rtmax = sqrt(safmax / 2);
            if (g1 > rtmin && g1 < rtmax) {
                //
                // Use unscaled algorithm
                //
                // The following two lines can be replaced by `d = abs( g )`.
                // This algorithm do not use the intrinsic complex abs.
                g2 = abssq(g);
                d = sqrt(g2);
                s = conj(g) / d;
                r = d;
            } else {
                //
                // Use scaled algorithm
                //
                u = min(safmax, max(safmin, g1));
                gs = g / u;
                // The following two lines can be replaced by `d = abs( gs )`.
                // This algorithm do not use the intrinsic complex abs.
                g2 = abssq(gs);
                d = sqrt(g2);
                s = conj(gs) / d;
                r = d * u;
            }
        }
    } else {
        f1 = max(abs(f.real()), abs(f.imag()));
        g1 = max(abs(g.real()), abs(g.imag()));
        rtmax = sqrt(safmax / 4);
        if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
            //
            // Use unscaled algorithm
            //
            f2 = abssq(f);
            g2 = abssq(g);
            h2 = f2 + g2;
            // safmin <= f2 <= h2 <= safmax
            if (f2 >= h2 * safmin) {
                // safmin <= f2/h2 <= 1, and h2/f2 is finite
                c = sqrt(f2 / h2);
                r = f / c;
                rtmax = rtmax * 2;
                if (f2 > rtmin && h2 < rtmax) {
                    // safmin <= sqrt( f2*h2 ) <= safmax
                    s = conj(g) * (f / sqrt(f2 * h2));
                } else {
                    s = conj(g) * (r / h2);
                }
            } else {
                // f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
                // Moreover,
                // safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
                // sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
                // Also,
                // g2 >> f2, which means that h2 = g2.
                d = sqrt(f2 * h2);
                c = f2 / d;
                if (c >= safmin) {
                    r = f / c;
                } else {
                    // f2 / sqrt(f2 * h2) < safmin, then
                    // sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
                    r = f * (h2 / d);
                }
                s = conj(g) * (f / d);
            }
        } else {
            //
            // Use scaled algorithm
            //
            u = min(safmax, max(safmin, f1, g1));
            gs = g / u;
            g2 = abssq(gs);
            if (f1 / u < rtmin) {
                //
                // f is not well-scaled when scaled by g1.
                // Use a different scaling for f.
                //
                v = min(safmax, max(safmin, f1));
                w = v / u;
                fs = f / v;
                f2 = abssq(fs);
                h2 = f2 * pow2(w) + g2;
            } else {
                //
                // Otherwise use the same scaling for f and g.
                //
                w = one;
                fs = f / u;
                f2 = abssq(fs);
                h2 = f2 + g2;
            }
            // safmin <= f2 <= h2 <= safmax
            if (f2 >= h2 * safmin) {
                // safmin <= f2/h2 <= 1, and h2/f2 is finite
                c = sqrt(f2 / h2);
                r = fs / c;
                rtmax = rtmax * 2;
                if (f2 > rtmin && h2 < rtmax) {
                    // safmin <= sqrt( f2*h2 ) <= safmax
                    s = conj(gs) * (fs / sqrt(f2 * h2));
                } else {
                    s = conj(gs) * (r / h2);
                }
            } else {
                // f2/h2 <= safmin may be subnormal, and h2/f2 may overflow.
                // Moreover,
                // safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= safmax,
                // sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax).
                // Also,
                // g2 >> f2, which means that h2 = g2.
                d = sqrt(f2 * h2);
                c = f2 / d;
                if (c >= safmin) {
                    r = fs / c;
                } else {
                    // f2 / sqrt(f2 * h2) < safmin, then
                    // sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2) <= h2 * (safmin / f2) <= h2 <= safmax
                    r = fs * (h2 / d);
                }
                s = conj(gs) * (fs / d);
            }
            // Rescale c and r
            c = c * w;
            r = r * u;
        }
    }
    a = r;
}
