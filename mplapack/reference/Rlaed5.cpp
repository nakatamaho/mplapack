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

void Rlaed5(INTEGER const i, REAL *d, REAL *z, REAL *delta, REAL const rho, REAL &dlam) {
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
    REAL del = d[2 - 1] - d[1 - 1];
    const REAL one = 1.0;
    const REAL two = 2.0;
    REAL w = 0.0;
    const REAL zero = 0.0;
    REAL b = 0.0;
    REAL c = 0.0;
    const REAL four = 4.0;
    REAL tau = 0.0;
    REAL temp = 0.0;
    if (i == 1) {
        w = one + two * rho * (z[2 - 1] * z[2 - 1] - z[1 - 1] * z[1 - 1]) / del;
        if (w > zero) {
            b = del + rho * (z[1 - 1] * z[1 - 1] + z[2 - 1] * z[2 - 1]);
            c = rho * z[1 - 1] * z[1 - 1] * del;
            //
            //           B > ZERO, always
            //
            tau = two * c / (b + sqrt(abs(b * b - four * c)));
            dlam = d[1 - 1] + tau;
            delta[1 - 1] = -z[1 - 1] / tau;
            delta[2 - 1] = z[2 - 1] / (del - tau);
        } else {
            b = -del + rho * (z[1 - 1] * z[1 - 1] + z[2 - 1] * z[2 - 1]);
            c = rho * z[2 - 1] * z[2 - 1] * del;
            if (b > zero) {
                tau = -two * c / (b + sqrt(b * b + four * c));
            } else {
                tau = (b - sqrt(b * b + four * c)) / two;
            }
            dlam = d[2 - 1] + tau;
            delta[1 - 1] = -z[1 - 1] / (del + tau);
            delta[2 - 1] = -z[2 - 1] / tau;
        }
        temp = sqrt(delta[1 - 1] * delta[1 - 1] + delta[2 - 1] * delta[2 - 1]);
        delta[1 - 1] = delta[1 - 1] / temp;
        delta[2 - 1] = delta[2 - 1] / temp;
    } else {
        //
        //     Now I=2
        //
        b = -del + rho * (z[1 - 1] * z[1 - 1] + z[2 - 1] * z[2 - 1]);
        c = rho * z[2 - 1] * z[2 - 1] * del;
        if (b > zero) {
            tau = (b + sqrt(b * b + four * c)) / two;
        } else {
            tau = two * c / (-b + sqrt(b * b + four * c));
        }
        dlam = d[2 - 1] + tau;
        delta[1 - 1] = -z[1 - 1] / (del + tau);
        delta[2 - 1] = -z[2 - 1] / tau;
        temp = sqrt(delta[1 - 1] * delta[1 - 1] + delta[2 - 1] * delta[2 - 1]);
        delta[1 - 1] = delta[1 - 1] / temp;
        delta[2 - 1] = delta[2 - 1] / temp;
    }
    //
    //     End OF Rlaed5
    //
}
