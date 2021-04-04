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

void Rlasq6(INTEGER const &i0, INTEGER const &n0, REAL *z, INTEGER const &pp, REAL &dmin, REAL &dmin1, REAL &dmin2, REAL &dn, REAL &dnm1, REAL &dnm2) {
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
    //     .. Parameter ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Function ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    if ((n0 - i0 - 1) <= 0) {
        return;
    }
    //
    REAL safmin = dlamch("Safe minimum");
    INTEGER j4 = 4 * i0 + pp - 3;
    REAL emin = z[(j4 + 4) - 1];
    REAL d = z[j4 - 1];
    dmin = d;
    //
    const REAL zero = 0.0;
    REAL temp = 0.0;
    if (pp == 0) {
        for (j4 = 4 * i0; j4 <= 4 * (n0 - 3); j4 = j4 + 4) {
            z[(j4 - 2) - 1] = d + z[(j4 - 1) - 1];
            if (z[(j4 - 2) - 1] == zero) {
                z[j4 - 1] = zero;
                d = z[(j4 + 1) - 1];
                dmin = d;
                emin = zero;
            } else if (safmin * z[(j4 + 1) - 1] < z[(j4 - 2) - 1] && safmin * z[(j4 - 2) - 1] < z[(j4 + 1) - 1]) {
                temp = z[(j4 + 1) - 1] / z[(j4 - 2) - 1];
                z[j4 - 1] = z[(j4 - 1) - 1] * temp;
                d = d * temp;
            } else {
                z[j4 - 1] = z[(j4 + 1) - 1] * (z[(j4 - 1) - 1] / z[(j4 - 2) - 1]);
                d = z[(j4 + 1) - 1] * (d / z[(j4 - 2) - 1]);
            }
            dmin = min(dmin, d);
            emin = min(emin, z[j4 - 1]);
        }
    } else {
        for (j4 = 4 * i0; j4 <= 4 * (n0 - 3); j4 = j4 + 4) {
            z[(j4 - 3) - 1] = d + z[j4 - 1];
            if (z[(j4 - 3) - 1] == zero) {
                z[(j4 - 1) - 1] = zero;
                d = z[(j4 + 2) - 1];
                dmin = d;
                emin = zero;
            } else if (safmin * z[(j4 + 2) - 1] < z[(j4 - 3) - 1] && safmin * z[(j4 - 3) - 1] < z[(j4 + 2) - 1]) {
                temp = z[(j4 + 2) - 1] / z[(j4 - 3) - 1];
                z[(j4 - 1) - 1] = z[j4 - 1] * temp;
                d = d * temp;
            } else {
                z[(j4 - 1) - 1] = z[(j4 + 2) - 1] * (z[j4 - 1] / z[(j4 - 3) - 1]);
                d = z[(j4 + 2) - 1] * (d / z[(j4 - 3) - 1]);
            }
            dmin = min(dmin, d);
            emin = min(emin, z[(j4 - 1) - 1]);
        }
    }
    //
    //     Unroll last two steps.
    //
    dnm2 = d;
    dmin2 = dmin;
    j4 = 4 * (n0 - 2) - pp;
    INTEGER j4p2 = j4 + 2 * pp - 1;
    z[(j4 - 2) - 1] = dnm2 + z[j4p2 - 1];
    if (z[(j4 - 2) - 1] == zero) {
        z[j4 - 1] = zero;
        dnm1 = z[(j4p2 + 2) - 1];
        dmin = dnm1;
        emin = zero;
    } else if (safmin * z[(j4p2 + 2) - 1] < z[(j4 - 2) - 1] && safmin * z[(j4 - 2) - 1] < z[(j4p2 + 2) - 1]) {
        temp = z[(j4p2 + 2) - 1] / z[(j4 - 2) - 1];
        z[j4 - 1] = z[j4p2 - 1] * temp;
        dnm1 = dnm2 * temp;
    } else {
        z[j4 - 1] = z[(j4p2 + 2) - 1] * (z[j4p2 - 1] / z[(j4 - 2) - 1]);
        dnm1 = z[(j4p2 + 2) - 1] * (dnm2 / z[(j4 - 2) - 1]);
    }
    dmin = min(dmin, dnm1);
    //
    dmin1 = dmin;
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    z[(j4 - 2) - 1] = dnm1 + z[j4p2 - 1];
    if (z[(j4 - 2) - 1] == zero) {
        z[j4 - 1] = zero;
        dn = z[(j4p2 + 2) - 1];
        dmin = dn;
        emin = zero;
    } else if (safmin * z[(j4p2 + 2) - 1] < z[(j4 - 2) - 1] && safmin * z[(j4 - 2) - 1] < z[(j4p2 + 2) - 1]) {
        temp = z[(j4p2 + 2) - 1] / z[(j4 - 2) - 1];
        z[j4 - 1] = z[j4p2 - 1] * temp;
        dn = dnm1 * temp;
    } else {
        z[j4 - 1] = z[(j4p2 + 2) - 1] * (z[j4p2 - 1] / z[(j4 - 2) - 1]);
        dn = z[(j4p2 + 2) - 1] * (dnm1 / z[(j4 - 2) - 1]);
    }
    dmin = min(dmin, dn);
    //
    z[(j4 + 2) - 1] = dn;
    z[(4 * n0 - pp) - 1] = emin;
    //
    //     End of Rlasq6
    //
}
