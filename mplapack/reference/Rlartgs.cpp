/*
 * Copyright (c) 2008-2021
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

#include <mpblas.h>
#include <mplapack.h>

void Rlartgs(REAL const x, REAL const y, REAL const sigma, REAL cs, REAL sn) {
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  ===================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     .. Executable Statements ..
    //
    REAL thresh = Rlamch("E");
    //
    //     Compute the first column of B**T*B - SIGMA^2*I, up to a scale
    //     factor.
    //
    const REAL zero = 0.0;
    REAL z = 0.0;
    REAL w = 0.0;
    const REAL one = 1.0;
    REAL s = 0.0;
    const REAL negone = -1.0;
    if ((sigma == zero && abs(x) < thresh) || (abs(x) == sigma && y == zero)) {
        z = zero;
        w = zero;
    } else if (sigma == zero) {
        if (x >= zero) {
            z = x;
            w = y;
        } else {
            z = -x;
            w = -y;
        }
    } else if (abs(x) < thresh) {
        z = -sigma * sigma;
        w = zero;
    } else {
        if (x >= zero) {
            s = one;
        } else {
            s = negone;
        }
        z = s * (abs(x) - sigma) * (s + sigma / x);
        w = s * y;
    }
    //
    //     Generate the rotation.
    //     CALL Rlartgp( Z, W, CS, SN, R ) might seem more natural;
    //     reordering the arguments ensures that if Z = 0 then the rotation
    //     is by PI/2.
    //
    REAL r = 0.0;
    Rlartgp(w, z, sn, cs, r);
    //
    //     End Rlartgs
    //
}
