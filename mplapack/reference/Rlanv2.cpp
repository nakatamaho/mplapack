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

void Rlanv2(REAL &a, REAL &b, REAL &c, REAL &d, REAL &rt1r, REAL &rt1i, REAL &rt2r, REAL &rt2i, REAL &cs, REAL &sn) {
    REAL safmin = 0.0;
    REAL eps = 0.0;
    const REAL two = 2.0;
    REAL safmn2 = 0.0;
    const REAL one = 1.0;
    REAL safmx2 = 0.0;
    const REAL zero = 0.0;
    REAL temp = 0.0;
    const REAL half = 0.5e+0;
    REAL p = 0.0;
    REAL bcmax = 0.0;
    REAL bcmis = 0.0;
    REAL scale = 0.0;
    REAL z = 0.0;
    const REAL multpl = 4.0e+0;
    REAL tau = 0.0;
    INTEGER count = 0;
    REAL sigma = 0.0;
    REAL aa = 0.0;
    REAL bb = 0.0;
    REAL cc = 0.0;
    REAL dd = 0.0;
    REAL sab = 0.0;
    REAL sac = 0.0;
    REAL cs1 = 0.0;
    REAL sn1 = 0.0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    safmin = Rlamch("S");
    eps = Rlamch("P");
    safmn2 = pow(Rlamch("B"), castINTEGER(log(safmin / eps) / log(Rlamch("B")) / two));
    safmx2 = one / safmn2;
    if (c == zero) {
        cs = one;
        sn = zero;
        //
    } else if (b == zero) {
        //
        //        Swap rows and columns
        //
        cs = zero;
        sn = one;
        temp = d;
        d = a;
        a = temp;
        b = -c;
        c = zero;
        //
    } else if ((a - d) == zero && sign(one, b) != sign(one, c)) {
        cs = one;
        sn = zero;
        //
    } else {
        //
        temp = a - d;
        p = half * temp;
        bcmax = max(abs(b), abs(c));
        bcmis = min(abs(b), abs(c)) * sign(one, b) * sign(one, c);
        scale = max(abs(p), bcmax);
        z = (p / scale) * p + (bcmax / scale) * bcmis;
        //
        //        If Z is of the order of the machine accuracy, postpone the
        //        decision on the nature of eigenvalues
        //
        if (z >= multpl * eps) {
            //
            //           Real eigenvalues. Compute A and D.
            //
            z = p + sign(sqrt(scale) * sqrt(z), p);
            a = d + z;
            d = d - (bcmax / z) * bcmis;
            //
            //           Compute B and the rotation matrix
            //
            tau = Rlapy2(c, z);
            cs = z / tau;
            sn = c / tau;
            b = b - c;
            c = zero;
            //
        } else {
            //
            //           Complex eigenvalues, or real (almost) equal eigenvalues.
            //           Make diagonal elements equal.
            //
            count = 0;
            sigma = b + c;
        statement_10:
            count++;
            scale = max(abs(temp), abs(sigma));
            if (scale >= safmx2) {
                sigma = sigma * safmn2;
                temp = temp * safmn2;
                if (count <= 20) {
                    goto statement_10;
                }
            }
            if (scale <= safmn2) {
                sigma = sigma * safmx2;
                temp = temp * safmx2;
                if (count <= 20) {
                    goto statement_10;
                }
            }
            p = half * temp;
            tau = Rlapy2(sigma, temp);
            cs = sqrt(half * (one + abs(sigma) / tau));
            sn = -(p / (tau * cs)) * sign(one, sigma);
            //
            //           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
            //                   [ CC  DD ]   [ C  D ] [ SN  CS ]
            //
            aa = a * cs + b * sn;
            bb = -a * sn + b * cs;
            cc = c * cs + d * sn;
            dd = -c * sn + d * cs;
            //
            //           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
            //                   [ C  D ]   [-SN  CS ] [ CC  DD ]
            //
            a = aa * cs + cc * sn;
            b = bb * cs + dd * sn;
            c = -aa * sn + cc * cs;
            d = -bb * sn + dd * cs;
            //
            temp = half * (a + d);
            a = temp;
            d = temp;
            //
            if (c != zero) {
                if (b != zero) {
                    if (sign(one, b) == sign(one, c)) {
                        //
                        //                    Real eigenvalues: reduce to upper triangular form
                        //
                        sab = sqrt(abs(b));
                        sac = sqrt(abs(c));
                        p = sign(sab * sac, c);
                        tau = one / sqrt(abs(b + c));
                        a = temp + p;
                        d = temp - p;
                        b = b - c;
                        c = zero;
                        cs1 = sab * tau;
                        sn1 = sac * tau;
                        temp = cs * cs1 - sn * sn1;
                        sn = cs * sn1 + sn * cs1;
                        cs = temp;
                    }
                } else {
                    b = -c;
                    c = zero;
                    temp = cs;
                    cs = -sn;
                    sn = temp;
                }
            }
        }
        //
    }
    //
    //     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
    //
    rt1r = a;
    rt2r = d;
    if (c == zero) {
        rt1i = zero;
        rt2i = zero;
    } else {
        rt1i = sqrt(abs(b)) * sqrt(abs(c));
        rt2i = -rt1i;
    }
    //
    //     End of Rlanv2
    //
}
