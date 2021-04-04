/*
 * Copyright (c) 2021
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

void Rlag2(REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL const &safmin, REAL &scale1, REAL &scale2, REAL &wr1, REAL &wr2, REAL &wi) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL rtmin = sqrt(safmin);
    const REAL one = 1.0;
    REAL rtmax = one / rtmin;
    REAL safmax = one / safmin;
    //
    //     Scale A
    //
    REAL anorm = max(abs(a[(1 - 1)]) + abs(a[(2 - 1)]), abs(a[(2 - 1) * lda]) + abs(a[(2 - 1) + (2 - 1) * lda]), safmin);
    REAL ascale = one / anorm;
    REAL a11 = ascale * a[(1 - 1)];
    REAL a21 = ascale * a[(2 - 1)];
    REAL a12 = ascale * a[(2 - 1) * lda];
    REAL a22 = ascale * a[(2 - 1) + (2 - 1) * lda];
    //
    //     Perturb B if necessary to insure non-singularity
    //
    REAL b11 = b[(1 - 1)];
    REAL b12 = b[(2 - 1) * ldb];
    REAL b22 = b[(2 - 1) + (2 - 1) * ldb];
    REAL bmin = rtmin * max(abs(b11), abs(b12), abs(b22), rtmin);
    if (abs(b11) < bmin) {
        b11 = sign[(bmin - 1) + (b11 - 1) * ldsign];
    }
    if (abs(b22) < bmin) {
        b22 = sign[(bmin - 1) + (b22 - 1) * ldsign];
    }
    //
    //     Scale B
    //
    REAL bnorm = max(abs(b11), abs(b12) + abs(b22), safmin);
    REAL bsize = max(abs(b11), abs(b22));
    REAL bscale = one / bsize;
    b11 = b11 * bscale;
    b12 = b12 * bscale;
    b22 = b22 * bscale;
    //
    //     Compute larger eigenvalue by method described by C. van Loan
    //
    //     ( AS is A shifted by -SHIFT*B )
    //
    REAL binv11 = one / b11;
    REAL binv22 = one / b22;
    REAL s1 = a11 * binv11;
    REAL s2 = a22 * binv22;
    REAL as12 = 0.0;
    REAL as22 = 0.0;
    REAL ss = 0.0;
    REAL abi22 = 0.0;
    const REAL two = 2.0e+0;
    const REAL half = one / two;
    REAL pp = 0.0;
    REAL shift = 0.0;
    REAL as11 = 0.0;
    if (abs(s1) <= abs(s2)) {
        as12 = a12 - s1 * b12;
        as22 = a22 - s1 * b22;
        ss = a21 * (binv11 * binv22);
        abi22 = as22 * binv22 - ss * b12;
        pp = half * abi22;
        shift = s1;
    } else {
        as12 = a12 - s2 * b12;
        as11 = a11 - s2 * b11;
        ss = a21 * (binv11 * binv22);
        abi22 = -ss * b12;
        pp = half * (as11 * binv11 + abi22);
        shift = s2;
    }
    REAL qq = ss * as12;
    REAL discr = 0.0;
    REAL r = 0.0;
    if (abs(pp * rtmin) >= one) {
        discr = pow2((rtmin * pp)) + qq * safmin;
        r = sqrt(abs(discr)) * rtmax;
    } else {
        if (pow2(pp) + abs(qq) <= safmin) {
            discr = pow2((rtmax * pp)) + qq * safmax;
            r = sqrt(abs(discr)) * rtmin;
        } else {
            discr = pow2(pp) + qq;
            r = sqrt(abs(discr));
        }
    }
    //
    //     Note: the test of R in the following IF is to cover the case when
    //           DISCR is small and negative and is flushed to zero during
    //           the calculation of R.  On machines which have a consistent
    //           flush-to-zero threshold and handle numbers above that
    //           threshold correctly, it would not be necessary.
    //
    const REAL zero = 0.0;
    REAL sum = 0.0;
    REAL diff = 0.0;
    REAL wbig = 0.0;
    REAL wsmall = 0.0;
    REAL wdet = 0.0;
    if (discr >= zero || r == zero) {
        sum = pp + sign[(r - 1) + (pp - 1) * ldsign];
        diff = pp - sign[(r - 1) + (pp - 1) * ldsign];
        wbig = shift + sum;
        //
        //        Compute smaller eigenvalue
        //
        wsmall = shift + diff;
        if (half * abs(wbig) > max(abs(wsmall), safmin)) {
            wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
            wsmall = wdet / wbig;
        }
        //
        //        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
        //        for WR1.
        //
        if (pp > abi22) {
            wr1 = min(wbig, wsmall);
            wr2 = max(wbig, wsmall);
        } else {
            wr1 = max(wbig, wsmall);
            wr2 = min(wbig, wsmall);
        }
        wi = zero;
    } else {
        //
        //        Complex eigenvalues
        //
        wr1 = shift + pp;
        wr2 = wr1;
        wi = r;
    }
    //
    //     Further scaling to avoid underflow and overflow in computing
    //     SCALE1 and overflow in computing w*B.
    //
    //     This scale factor (WSCALE) is bounded from above using C1 and C2,
    //     and from below using C3 and C4.
    //        C1 implements the condition  s A  must never overflow.
    //        C2 implements the condition  w B  must never overflow.
    //        C3, with C2,
    //           implement the condition that s A - w B must never overflow.
    //        C4 implements the condition  s    should not underflow.
    //        C5 implements the condition  max(s,|w|) should be at least 2.
    //
    REAL c1 = bsize * (safmin * max(one, ascale));
    REAL c2 = safmin * max(one, bnorm);
    REAL c3 = bsize * safmin;
    REAL c4 = 0.0;
    if (ascale <= one && bsize <= one) {
        c4 = min(one, (ascale / safmin) * bsize);
    } else {
        c4 = one;
    }
    REAL c5 = 0.0;
    if (ascale <= one || bsize <= one) {
        c5 = min(one, ascale * bsize);
    } else {
        c5 = one;
    }
    //
    //     Scale first eigenvalue
    //
    REAL wabs = abs(wr1) + abs(wi);
    const REAL fuzzy1 = one + 1.0e-5;
    REAL wsize = max(safmin, c1, fuzzy1 * (wabs * c2 + c3), min(c4, half * max(wabs, c5)));
    REAL wscale = 0.0;
    if (wsize != one) {
        wscale = one / wsize;
        if (wsize > one) {
            scale1 = (max(ascale, bsize) * wscale) * min(ascale, bsize);
        } else {
            scale1 = (min(ascale, bsize) * wscale) * max(ascale, bsize);
        }
        wr1 = wr1 * wscale;
        if (wi != zero) {
            wi = wi * wscale;
            wr2 = wr1;
            scale2 = scale1;
        }
    } else {
        scale1 = ascale * bsize;
        scale2 = scale1;
    }
    //
    //     Scale second eigenvalue (if real)
    //
    if (wi == zero) {
        wsize = max(safmin, c1, fuzzy1 * (abs(wr2) * c2 + c3), min(c4, half * max(abs(wr2), c5)));
        if (wsize != one) {
            wscale = one / wsize;
            if (wsize > one) {
                scale2 = (max(ascale, bsize) * wscale) * min(ascale, bsize);
            } else {
                scale2 = (min(ascale, bsize) * wscale) * max(ascale, bsize);
            }
            wr2 = wr2 * wscale;
        } else {
            scale2 = ascale * bsize;
        }
    }
    //
    //     End of Rlag2
    //
}
