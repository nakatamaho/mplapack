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

#include <mpblas.h>

struct common : fem::common {
    fem::cmn_sve drotmg_sve;

    common(int argc, char const *argv[]) : fem::common(argc, argv) {}
};

struct drotmg_save {
    double gam;
    double gamsq;
    double one;
    double rgamsq;
    double two;
    double zero;

    drotmg_save() : gam(0.0), gamsq(0.0), one(0.0), rgamsq(0.0), two(0.0), zero(0.0) {}
};

// > \brief \b DROTMG
// >
// > \verbatim
// >
// >    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
// >    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T.
// >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
// >
// >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
// >
// >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
// >    H=(          )    (          )    (          )    (          )
// >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
// >    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
// >    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
// >    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
// >
// >    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
// >    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
// >    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
// >
// > \endverbatim
//
// Arguments:
// > \param[in,out] DD1
// > \verbatim
// >          DD1 is DOUBLE PRECISION
// > \endverbatim
// >
// > \param[in,out] DD2
// > \verbatim
// >          DD2 is DOUBLE PRECISION
// > \endverbatim
// >
// > \param[in,out] DX1
// > \verbatim
// >          DX1 is DOUBLE PRECISION
// > \endverbatim
// >
// > \param[in] DY1
// > \verbatim
// >          DY1 is DOUBLE PRECISION
// > \endverbatim
// >
// > \param[out] DPARAM
// > \verbatim
// >          DPARAM is DOUBLE PRECISION array, dimension (5)
// >     DPARAM(1)=DFLAG
// >     DPARAM(2)=DH11
// >     DPARAM(3)=DH21
// >     DPARAM(4)=DH12
// >     DPARAM(5)=DH22
// > \endverbatim
//
// Authors:
void Rrotmg(common &cmn, REAL const &dd1, REAL const &dd2, REAL const &dx1, REAL const &dy1, REAL *dparam) {
    FEM_CMN_SVE(drotmg);
    // SAVE
    double &gam = sve.gam;
    double &gamsq = sve.gamsq;
    double &one = sve.one;
    double &rgamsq = sve.rgamsq;
    double &two = sve.two;
    double &zero = sve.zero;
    //
    if (is_called_first_time) {
        zero = 0.0;
        one = 1.0;
        two = 2.0;
        gam = 4096.0;
        gamsq = 16777216.0;
        rgamsq = 5.960464500000001e-08;
    }
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
        // GO ZERO-H-D-AND-DX1..
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
        // CASE-DD1-NONNEGATIVE
        dp2 = dd2 * dy1;
        if (dp2 == zero) {
            dflag = -two;
            dparam[1 - 1] = dflag;
            return;
        }
        // REGULAR-CASE..
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
                // This code path if here for safety. We do not expect this
                // condition to ever hold except in edge cases with rounding
                // errors. See DOI: 10.1145/355841.355847
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
                // GO ZERO-H-D-AND-DX1..
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
        // PROCEDURE..SCALE-CHECK
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
