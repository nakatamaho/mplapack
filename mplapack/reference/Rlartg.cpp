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

// Derived from LAPACK routine DLARTG.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_arithmetic_params.h>

void Rlartg(REAL const f, REAL const g, REAL &c, REAL &s, REAL &r) {
    const auto &ap = mplapack::get_arithmetic_params<REAL>();
    const REAL safmin = ap.safmin;
    REAL rtmin = sqrt(safmin);
    const REAL safmax = ap.safmax;
    REAL rtmax = sqrt(safmax / 2);
    //
    REAL f1 = abs(f);
    REAL g1 = abs(g);
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL d = 0.0;
    REAL u = 0.0;
    REAL fs = 0.0;
    REAL gs = 0.0;
    if (g == zero) {
        c = one;
        s = zero;
        r = f;
    } else if (f == zero) {
        c = zero;
        s = sign(one, g);
        r = g1;
    } else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
        d = sqrt(f * f + g * g);
        c = f1 / d;
        r = sign(d, f);
        s = g / r;
    } else {
        u = min(safmax, max(safmin, f1, g1));
        fs = f / u;
        gs = g / u;
        d = sqrt(fs * fs + gs * gs);
        c = abs(fs) / d;
        r = sign(d, f);
        s = gs / r;
        r = r * u;
    }
}
