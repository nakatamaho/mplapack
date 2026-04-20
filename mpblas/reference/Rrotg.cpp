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

// Derived from BLAS routine DROTG.
// Original BLAS authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack_arithmetic_params.h>

void Rrotg(REAL &a, REAL &b, REAL &c, REAL &s) {
    //
    REAL anorm = abs(a);
    REAL bnorm = abs(b);
    const REAL zero = 0.0;
    const REAL one = 1.0;
    const auto &ap = mplapack::get_arithmetic_params<REAL>();
    const REAL safmin = ap.safmin;
    const REAL safmax = ap.safmax;
    REAL scl = 0.0;
    REAL sigma = 0.0;
    REAL r = 0.0;
    REAL z = 0.0;
    if (bnorm == zero) {
        c = one;
        s = zero;
        b = zero;
    } else if (anorm == zero) {
        c = zero;
        s = one;
        a = b;
        b = one;
    } else {
        scl = min(safmax, max(safmin, anorm, bnorm));
        if (anorm > bnorm) {
            sigma = sign(one, a);
        } else {
            sigma = sign(one, b);
        }
        r = sigma * (scl * sqrt(pow2((a / scl)) + pow2((b / scl))));
        c = a / r;
        s = b / r;
        if (anorm > bnorm) {
            z = s;
        } else if (c != zero) {
            z = one / c;
        } else {
            z = one;
        }
        a = r;
        b = z;
    }
}
