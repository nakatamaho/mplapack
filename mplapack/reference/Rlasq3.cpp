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

void Rlasq3(INTEGER const i0, INTEGER &n0, REAL *z, INTEGER &pp, REAL &dmin, REAL &sigma, REAL &desig, REAL &qmax, INTEGER &nfail, INTEGER &iter, INTEGER &ndiv, bool const ieee, INTEGER &ttype, REAL dmin1, REAL &dmin2, REAL dn, REAL dn1, REAL dn2, REAL g, REAL &tau) {
    INTEGER n0in = 0;
    REAL eps = 0.0;
    const REAL hundrd = 100.0;
    REAL tol = 0.0;
    REAL tol2 = 0.0;
    INTEGER nn = 0;
    REAL s = 0.0;
    const REAL half = 0.5e0;
    REAL t = 0.0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    const REAL cbias = 1.50e0;
    INTEGER ipn4 = 0;
    INTEGER j4 = 0;
    REAL temp = 0.0;
    const REAL two = 2.0;
    const REAL qurtr = 0.250e0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Function ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    n0in = n0;
    eps = Rlamch("Precision");
    tol = eps * hundrd;
    tol2 = pow2(tol);
//
//     Check for deflation.
//
statement_10:
    //
    if (n0 < i0) {
        return;
    }
    if (n0 == i0) {
        goto statement_20;
    }
    nn = 4 * n0 + pp;
    if (n0 == (i0 + 1)) {
        goto statement_40;
    }
    //
    //     Check whether E(N0-1) is negligible, 1 eigenvalue.
    //
    if (z[(nn - 5) - 1] > tol2 * (sigma + z[(nn - 3) - 1]) && z[(nn - 2 * pp - 4) - 1] > tol2 * z[(nn - 7) - 1]) {
        goto statement_30;
    }
//
statement_20:
    //
    z[(4 * n0 - 3) - 1] = z[(4 * n0 + pp - 3) - 1] + sigma;
    n0 = n0 - 1;
    goto statement_10;
//
//     Check  whether E(N0-2) is negligible, 2 eigenvalues.
//
statement_30:
    //
    if (z[(nn - 9) - 1] > tol2 * sigma && z[(nn - 2 * pp - 8) - 1] > tol2 * z[(nn - 11) - 1]) {
        goto statement_50;
    }
//
statement_40:
    //
    if (z[(nn - 3) - 1] > z[(nn - 7) - 1]) {
        s = z[(nn - 3) - 1];
        z[(nn - 3) - 1] = z[(nn - 7) - 1];
        z[(nn - 7) - 1] = s;
    }
    t = half * ((z[(nn - 7) - 1] - z[(nn - 3) - 1]) + z[(nn - 5) - 1]);
    if (z[(nn - 5) - 1] > z[(nn - 3) - 1] * tol2 && t != zero) {
        s = z[(nn - 3) - 1] * (z[(nn - 5) - 1] / t);
        if (s <= t) {
            s = z[(nn - 3) - 1] * (z[(nn - 5) - 1] / (t * (one + sqrt(one + s / t))));
        } else {
            s = z[(nn - 3) - 1] * (z[(nn - 5) - 1] / (t + sqrt(t) * sqrt(t + s)));
        }
        t = z[(nn - 7) - 1] + (s + z[(nn - 5) - 1]);
        z[(nn - 3) - 1] = z[(nn - 3) - 1] * (z[(nn - 7) - 1] / t);
        z[(nn - 7) - 1] = t;
    }
    z[(4 * n0 - 7) - 1] = z[(nn - 7) - 1] + sigma;
    z[(4 * n0 - 3) - 1] = z[(nn - 3) - 1] + sigma;
    n0 = n0 - 2;
    goto statement_10;
//
statement_50:
    if (pp == 2) {
        pp = 0;
    }
    //
    //     Reverse the qd-array, if warranted.
    //
    if (dmin <= zero || n0 < n0in) {
        if (cbias * z[(4 * i0 + pp - 3) - 1] < z[(4 * n0 + pp - 3) - 1]) {
            ipn4 = 4 * (i0 + n0);
            for (j4 = 4 * i0; j4 <= 2 * (i0 + n0 - 1); j4 = j4 + 4) {
                temp = z[(j4 - 3) - 1];
                z[(j4 - 3) - 1] = z[(ipn4 - j4 - 3) - 1];
                z[(ipn4 - j4 - 3) - 1] = temp;
                temp = z[(j4 - 2) - 1];
                z[(j4 - 2) - 1] = z[(ipn4 - j4 - 2) - 1];
                z[(ipn4 - j4 - 2) - 1] = temp;
                temp = z[(j4 - 1) - 1];
                z[(j4 - 1) - 1] = z[(ipn4 - j4 - 5) - 1];
                z[(ipn4 - j4 - 5) - 1] = temp;
                temp = z[j4 - 1];
                z[j4 - 1] = z[(ipn4 - j4 - 4) - 1];
                z[(ipn4 - j4 - 4) - 1] = temp;
            }
            if (n0 - i0 <= 4) {
                z[(4 * n0 + pp - 1) - 1] = z[(4 * i0 + pp - 1) - 1];
                z[(4 * n0 - pp) - 1] = z[(4 * i0 - pp) - 1];
            }
            dmin2 = min(dmin2, z[(4 * n0 + pp - 1) - 1]);
            z[(4 * n0 + pp - 1) - 1] = min({z[(4 * n0 + pp - 1) - 1], z[(4 * i0 + pp - 1) - 1], z[(4 * i0 + pp + 3) - 1]});
            z[(4 * n0 - pp) - 1] = min({z[(4 * n0 - pp) - 1], z[(4 * i0 - pp) - 1], z[(4 * i0 - pp + 4) - 1]});
            qmax = max({qmax, z[(4 * i0 + pp - 3) - 1], z[(4 * i0 + pp + 1) - 1]});
            dmin = -zero;
        }
    }
    //
    //     Choose a shift.
    //
    Rlasq4(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
//
//     Call dqds until DMIN > 0.
//
statement_70:
    //
    Rlasq5(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dn1, dn2, ieee, eps);
    //
    ndiv += (n0 - i0 + 2);
    iter++;
    //
    //     Check status.
    //
    if (dmin >= zero && dmin1 >= zero) {
        //
        //        Success.
        //
        goto statement_90;
        //
    } else if (dmin < zero && dmin1 > zero && z[(4 * (n0 - 1) - pp) - 1] < tol * (sigma + dn1) && abs(dn) < tol * sigma) {
        //
        //        Convergence hidden by negative DN.
        //
        z[(4 * (n0 - 1) - pp + 2) - 1] = zero;
        dmin = zero;
        goto statement_90;
    } else if (dmin < zero) {
        //
        //        TAU too big. Select new TAU and try again.
        //
        nfail++;
        if (ttype < -22) {
            //
            //           Failed twice. Play it safe.
            //
            tau = zero;
        } else if (dmin1 > zero) {
            //
            //           Late failure. Gives excellent shift.
            //
            tau = (tau + dmin) * (one - two * eps);
            ttype = ttype - 11;
        } else {
            //
            //           Early failure. Divide by 4.
            //
            tau = qurtr * tau;
            ttype = ttype - 12;
        }
        goto statement_70;
    } else if (Risnan(dmin)) {
        //
        //        NaN.
        //
        if (tau == zero) {
            goto statement_80;
        } else {
            tau = zero;
            goto statement_70;
        }
    } else {
        //
        //        Possible underflow. Play it safe.
        //
        goto statement_80;
    }
//
//     Risk of underflow.
//
statement_80:
    Rlasq6(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dn1, dn2);
    ndiv += (n0 - i0 + 2);
    iter++;
    tau = zero;
//
statement_90:
    if (tau < sigma) {
        desig += tau;
        t = sigma + desig;
        desig = desig - (t - sigma);
    } else {
        t = sigma + tau;
        desig += sigma - (t - tau);
    }
    sigma = t;
    //
    //     End of Rlasq3
    //
}
