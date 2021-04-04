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

void Rlasq4(INTEGER const &i0, INTEGER const &n0, REAL *z, INTEGER const &pp, INTEGER const &n0in, REAL const &dmin, REAL const &dmin1, REAL const &dmin2, REAL const &dn, REAL const &dn1, REAL const &dn2, REAL &tau, INTEGER &ttype, REAL &g) {
    const REAL zero = 0.0;
    INTEGER nn = 0;
    REAL b1 = 0.0;
    REAL b2 = 0.0;
    REAL a2 = 0.0;
    const REAL qurtr = 0.250e0;
    REAL gap2 = 0.0;
    REAL gap1 = 0.0;
    const REAL half = 0.50e0;
    REAL s = 0.0;
    const REAL third = 0.3330e0;
    REAL gam = 0.0;
    INTEGER np = 0;
    INTEGER i4 = 0;
    const REAL hundrd = 100.0;
    const REAL cnst1 = 0.5630e0;
    const REAL cnst3 = 1.050e0;
    const REAL one = 1.0;
    const REAL cnst2 = 1.010e0;
    const REAL two = 2.0;
    //
    //  -- LAPACK computational routine --
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
    //     A negative DMIN forces the shift to take that absolute value
    //     TTYPE records the type of shift.
    //
    if (dmin <= zero) {
        tau = -dmin;
        ttype = -1;
        return;
    }
    //
    nn = 4 * n0 + pp;
    if (n0in == n0) {
        //
        //        No eigenvalues deflated.
        //
        if (dmin == dn || dmin == dn1) {
            //
            b1 = sqrt(z[(nn - 3) - 1]) * sqrt(z[(nn - 5) - 1]);
            b2 = sqrt(z[(nn - 7) - 1]) * sqrt(z[(nn - 9) - 1]);
            a2 = z[(nn - 7) - 1] + z[(nn - 5) - 1];
            //
            //           Cases 2 and 3.
            //
            if (dmin == dn && dmin1 == dn1) {
                gap2 = dmin2 - a2 - dmin2 * qurtr;
                if (gap2 > zero && gap2 > b2) {
                    gap1 = a2 - dn - (b2 / gap2) * b2;
                } else {
                    gap1 = a2 - dn - (b1 + b2);
                }
                if (gap1 > zero && gap1 > b1) {
                    s = max(dn - (b1 / gap1) * b1, half * dmin);
                    ttype = -2;
                } else {
                    s = zero;
                    if (dn > b1) {
                        s = dn - b1;
                    }
                    if (a2 > (b1 + b2)) {
                        s = min(s, a2 - (b1 + b2));
                    }
                    s = max(s, third * dmin);
                    ttype = -3;
                }
            } else {
                //
                //              Case 4.
                //
                ttype = -4;
                s = qurtr * dmin;
                if (dmin == dn) {
                    gam = dn;
                    a2 = zero;
                    if (z[(nn - 5) - 1] > z[(nn - 7) - 1]) {
                        return;
                    }
                    b2 = z[(nn - 5) - 1] / z[(nn - 7) - 1];
                    np = nn - 9;
                } else {
                    np = nn - 2 * pp;
                    gam = dn1;
                    if (z[(np - 4) - 1] > z[(np - 2) - 1]) {
                        return;
                    }
                    a2 = z[(np - 4) - 1] / z[(np - 2) - 1];
                    if (z[(nn - 9) - 1] > z[(nn - 11) - 1]) {
                        return;
                    }
                    b2 = z[(nn - 9) - 1] / z[(nn - 11) - 1];
                    np = nn - 13;
                }
                //
                //              Approximate contribution to norm squared from I < NN-1.
                //
                a2 += b2;
                for (i4 = np; i4 >= 4 * i0 - 1 + pp; i4 = i4 - 4) {
                    if (b2 == zero) {
                        goto statement_20;
                    }
                    b1 = b2;
                    if (z[i4 - 1] > z[(i4 - 2) - 1]) {
                        return;
                    }
                    b2 = b2 * (z[i4 - 1] / z[(i4 - 2) - 1]);
                    a2 += b2;
                    if (hundrd * max(b2, b1) < a2 || cnst1 < a2) {
                        goto statement_20;
                    }
                }
            statement_20:
                a2 = cnst3 * a2;
                //
                //              Rayleigh quotient residual bound.
                //
                if (a2 < cnst1) {
                    s = gam * (one - sqrt(a2)) / (one + a2);
                }
            }
        } else if (dmin == dn2) {
            //
            //           Case 5.
            //
            ttype = -5;
            s = qurtr * dmin;
            //
            //           Compute contribution to norm squared from I > NN-2.
            //
            np = nn - 2 * pp;
            b1 = z[(np - 2) - 1];
            b2 = z[(np - 6) - 1];
            gam = dn2;
            if (z[(np - 8) - 1] > b2 || z[(np - 4) - 1] > b1) {
                return;
            }
            a2 = (z[(np - 8) - 1] / b2) * (one + z[(np - 4) - 1] / b1);
            //
            //           Approximate contribution to norm squared from I < NN-2.
            //
            if (n0 - i0 > 2) {
                b2 = z[(nn - 13) - 1] / z[(nn - 15) - 1];
                a2 += b2;
                for (i4 = nn - 17; i4 >= 4 * i0 - 1 + pp; i4 = i4 - 4) {
                    if (b2 == zero) {
                        goto statement_40;
                    }
                    b1 = b2;
                    if (z[i4 - 1] > z[(i4 - 2) - 1]) {
                        return;
                    }
                    b2 = b2 * (z[i4 - 1] / z[(i4 - 2) - 1]);
                    a2 += b2;
                    if (hundrd * max(b2, b1) < a2 || cnst1 < a2) {
                        goto statement_40;
                    }
                }
            statement_40:
                a2 = cnst3 * a2;
            }
            //
            if (a2 < cnst1) {
                s = gam * (one - sqrt(a2)) / (one + a2);
            }
        } else {
            //
            //           Case 6, no information to guide us.
            //
            if (ttype == -6) {
                g += third * (one - g);
            } else if (ttype == -18) {
                g = qurtr * third;
            } else {
                g = qurtr;
            }
            s = g * dmin;
            ttype = -6;
        }
        //
    } else if (n0in == (n0 + 1)) {
        //
        //        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
        //
        if (dmin1 == dn1 && dmin2 == dn2) {
            //
            //           Cases 7 and 8.
            //
            ttype = -7;
            s = third * dmin1;
            if (z[(nn - 5) - 1] > z[(nn - 7) - 1]) {
                return;
            }
            b1 = z[(nn - 5) - 1] / z[(nn - 7) - 1];
            b2 = b1;
            if (b2 == zero) {
                goto statement_60;
            }
            for (i4 = 4 * n0 - 9 + pp; i4 >= 4 * i0 - 1 + pp; i4 = i4 - 4) {
                a2 = b1;
                if (z[i4 - 1] > z[(i4 - 2) - 1]) {
                    return;
                }
                b1 = b1 * (z[i4 - 1] / z[(i4 - 2) - 1]);
                b2 += b1;
                if (hundrd * max(b1, a2) < b2) {
                    goto statement_60;
                }
            }
        statement_60:
            b2 = sqrt(cnst3 * b2);
            a2 = dmin1 / (one + pow2(b2));
            gap2 = half * dmin2 - a2;
            if (gap2 > zero && gap2 > b2 * a2) {
                s = max(s, a2 * (one - cnst2 * a2 * (b2 / gap2) * b2));
            } else {
                s = max(s, a2 * (one - cnst2 * b2));
                ttype = -8;
            }
        } else {
            //
            //           Case 9.
            //
            s = qurtr * dmin1;
            if (dmin1 == dn1) {
                s = half * dmin1;
            }
            ttype = -9;
        }
        //
    } else if (n0in == (n0 + 2)) {
        //
        //        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
        //
        //        Cases 10 and 11.
        //
        if (dmin2 == dn2 && two * z[(nn - 5) - 1] < z[(nn - 7) - 1]) {
            ttype = -10;
            s = third * dmin2;
            if (z[(nn - 5) - 1] > z[(nn - 7) - 1]) {
                return;
            }
            b1 = z[(nn - 5) - 1] / z[(nn - 7) - 1];
            b2 = b1;
            if (b2 == zero) {
                goto statement_80;
            }
            for (i4 = 4 * n0 - 9 + pp; i4 >= 4 * i0 - 1 + pp; i4 = i4 - 4) {
                if (z[i4 - 1] > z[(i4 - 2) - 1]) {
                    return;
                }
                b1 = b1 * (z[i4 - 1] / z[(i4 - 2) - 1]);
                b2 += b1;
                if (hundrd * b1 < b2) {
                    goto statement_80;
                }
            }
        statement_80:
            b2 = sqrt(cnst3 * b2);
            a2 = dmin2 / (one + pow2(b2));
            gap2 = z[(nn - 7) - 1] + z[(nn - 9) - 1] - sqrt(z[(nn - 11) - 1]) * sqrt(z[(nn - 9) - 1]) - a2;
            if (gap2 > zero && gap2 > b2 * a2) {
                s = max(s, a2 * (one - cnst2 * a2 * (b2 / gap2) * b2));
            } else {
                s = max(s, a2 * (one - cnst2 * b2));
            }
        } else {
            s = qurtr * dmin2;
            ttype = -11;
        }
    } else if (n0in > (n0 + 2)) {
        //
        //        Case 12, more than two eigenvalues deflated. No information.
        //
        s = zero;
        ttype = -12;
    }
    //
    tau = s;
    //
    //     End of Rlasq4
    //
}
