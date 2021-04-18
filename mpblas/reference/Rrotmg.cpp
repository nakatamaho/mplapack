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

void Rrotmg(common &cmn, REAL &dd1, REAL &dd2, REAL &dx1, REAL const dy1, REAL *dparam) {
    FEM_CMN_SVE(Rrotmg);
    // SAVE
    REAL &gam = sve.gam;
    REAL &gamsq = sve.gamsq;
    REAL &one = sve.one;
    REAL &rgamsq = sve.rgamsq;
    REAL &two = sve.two;
    REAL &zero = sve.zero;
    //
    if (is_called_first_time) {
        zero = 0.0;
        one = 1.0;
        two = 2.0;
        gam = 4096.0;
        gamsq = 16777216;
        rgamsq = 5.9604645e-8;
    }
    //
    //  -- Reference BLAS level1 routine --
    //  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //
    //     ..
    //
    REAL dflag = 0.0;
    REAL dh11 = 0.0;
    REAL dh12 = 0.0;
    REAL dh21 = 0.0;
    REAL dh22 = 0.0;
    REAL dp2 = 0.0;
    REAL dp1 = 0.0;
    REAL dq2 = 0.0;
    REAL dq1 = 0.0;
    REAL du = 0.0;
    REAL dtemp = 0.0;
    if (dd1 < zero) {
        //        GO ZERO-H-D-AND-DX1..
        dflag = -one;
        dh11 = zero;
        dh12 = zero;
        dh21 = zero;
        dh22 = zero;
        //
        dd1 = zero;
        dd2 = zero;
        dx1 = zero;
    } else {
        //        CASE-DD1-NONNEGATIVE
        dp2 = dd2 * dy1;
        if (dp2 == zero) {
            dflag = -two;
            dparam[1 - 1] = dflag;
            return;
        }
        //        REGULAR-CASE..
        dp1 = dd1 * dx1;
        dq2 = dp2 * dy1;
        dq1 = dp1 * dx1;
        //
        if (abs(dq1) > abs(dq2)) {
            dh21 = -dy1 / dx1;
            dh12 = dp2 / dp1;
            //
            du = one - dh12 * dh21;
            //
            if (du > zero) {
                dflag = zero;
                dd1 = dd1 / du;
                dd2 = dd2 / du;
                dx1 = dx1 * du;
            } else {
                //            This code path if here for safety. We do not expect this
                //            condition to ever hold except in edge cases with rounding
                //            errors. See DOI: 10.1145/355841.355847
                dflag = -one;
                dh11 = zero;
                dh12 = zero;
                dh21 = zero;
                dh22 = zero;
                //
                dd1 = zero;
                dd2 = zero;
                dx1 = zero;
            }
        } else {
            //
            if (dq2 < zero) {
                //              GO ZERO-H-D-AND-DX1..
                dflag = -one;
                dh11 = zero;
                dh12 = zero;
                dh21 = zero;
                dh22 = zero;
                //
                dd1 = zero;
                dd2 = zero;
                dx1 = zero;
            } else {
                dflag = one;
                dh11 = dp1 / dp2;
                dh22 = dx1 / dy1;
                du = one + dh11 * dh22;
                dtemp = dd2 / du;
                dd2 = dd1 / du;
                dd1 = dtemp;
                dx1 = dy1 * du;
            }
        }
        //
        //     PROCEDURE..SCALE-CHECK
        if (dd1 != zero) {
            while ((dd1 <= rgamsq) || (dd1 >= gamsq)) {
                if (dflag == zero) {
                    dh11 = one;
                    dh22 = one;
                    dflag = -one;
                } else {
                    dh21 = -one;
                    dh12 = one;
                    dflag = -one;
                }
                if (dd1 <= rgamsq) {
                    dd1 = dd1 * pow2(gam);
                    dx1 = dx1 / gam;
                    dh11 = dh11 / gam;
                    dh12 = dh12 / gam;
                } else {
                    dd1 = dd1 / pow2(gam);
                    dx1 = dx1 * gam;
                    dh11 = dh11 * gam;
                    dh12 = dh12 * gam;
                }
            }
        }
        //
        if (dd2 != zero) {
            while ((abs(dd2) <= rgamsq) || (abs(dd2) >= gamsq)) {
                if (dflag == zero) {
                    dh11 = one;
                    dh22 = one;
                    dflag = -one;
                } else {
                    dh21 = -one;
                    dh12 = one;
                    dflag = -one;
                }
                if (abs(dd2) <= rgamsq) {
                    dd2 = dd2 * pow2(gam);
                    dh21 = dh21 / gam;
                    dh22 = dh22 / gam;
                } else {
                    dd2 = dd2 / pow2(gam);
                    dh21 = dh21 * gam;
                    dh22 = dh22 * gam;
                }
            }
        }
        //
    }
    //
    if (dflag < zero) {
        dparam[2 - 1] = dh11;
        dparam[3 - 1] = dh21;
        dparam[4 - 1] = dh12;
        dparam[5 - 1] = dh22;
    } else if (dflag == zero) {
        dparam[3 - 1] = dh21;
        dparam[4 - 1] = dh12;
    } else {
        dparam[2 - 1] = dh11;
        dparam[5 - 1] = dh22;
    }
    //
    dparam[1 - 1] = dflag;
}
