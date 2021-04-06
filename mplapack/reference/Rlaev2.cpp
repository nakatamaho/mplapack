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

void Rlaev2(REAL const a, REAL const b, REAL const c, REAL &rt1, REAL &rt2, REAL &cs1, REAL &sn1) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    // =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Compute the eigenvalues
    //
    REAL sm = a + c;
    REAL df = a - c;
    REAL adf = abs(df);
    REAL tb = b + b;
    REAL ab = abs(tb);
    REAL acmx = 0.0;
    REAL acmn = 0.0;
    if (abs(a) > abs(c)) {
        acmx = a;
        acmn = c;
    } else {
        acmx = c;
        acmn = a;
    }
    const REAL one = 1.0;
    REAL rt = 0.0;
    const REAL two = 2.0;
    if (adf > ab) {
        rt = adf * sqrt(one + pow2((ab / adf)));
    } else if (adf < ab) {
        rt = ab * sqrt(one + pow2((adf / ab)));
    } else {
        //
        //        Includes case AB=ADF=0
        //
        rt = ab * sqrt(two);
    }
    const REAL zero = 0.0;
    const REAL half = 0.5e0;
    INTEGER sgn1 = 0;
    if (sm < zero) {
        rt1 = half * (sm - rt);
        sgn1 = -1;
        //
        //        Order of execution important.
        //        To get fully accurate smaller eigenvalue,
        //        next line needs to be executed in higher precision.
        //
        rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
    } else if (sm > zero) {
        rt1 = half * (sm + rt);
        sgn1 = 1;
        //
        //        Order of execution important.
        //        To get fully accurate smaller eigenvalue,
        //        next line needs to be executed in higher precision.
        //
        rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
    } else {
        //
        //        Includes case RT1 = RT2 = 0
        //
        rt1 = half * rt;
        rt2 = -half * rt;
        sgn1 = 1;
    }
    //
    //     Compute the eigenvector
    //
    REAL cs = 0.0;
    INTEGER sgn2 = 0;
    if (df >= zero) {
        cs = df + rt;
        sgn2 = 1;
    } else {
        cs = df - rt;
        sgn2 = -1;
    }
    REAL acs = abs(cs);
    REAL ct = 0.0;
    REAL tn = 0.0;
    if (acs > ab) {
        ct = -tb / cs;
        sn1 = one / sqrt(one + ct * ct);
        cs1 = ct * sn1;
    } else {
        if (ab == zero) {
            cs1 = one;
            sn1 = zero;
        } else {
            tn = -cs / tb;
            cs1 = one / sqrt(one + tn * tn);
            sn1 = tn * cs1;
        }
    }
    if (sgn1 == sgn2) {
        tn = cs1;
        cs1 = -sn1;
        sn1 = tn;
    }
    //
    //     End of Rlaev2
    //
}
