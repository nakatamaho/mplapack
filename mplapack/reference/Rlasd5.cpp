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

void Rlasd5(INTEGER const i, REAL *d, REAL *z, REAL *delta, REAL const rho, REAL &dsigma, REAL *work) {
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
    REAL del = d[2 - 1] - d[1 - 1];
    REAL delsq = del * (d[2 - 1] + d[1 - 1]);
    const REAL one = 1.0;
    const REAL four = 4.0e+0;
    const REAL three = 3.0e+0;
    REAL w = 0.0;
    const REAL zero = 0.0;
    REAL b = 0.0;
    REAL c = 0.0;
    const REAL two = 2.0e+0;
    REAL tau = 0.0;
    if (i == 1) {
        w = one + four * rho * (z[2 - 1] * z[2 - 1] / (d[1 - 1] + three * d[2 - 1]) - z[1 - 1] * z[1 - 1] / (three * d[1 - 1] + d[2 - 1])) / del;
        if (w > zero) {
            b = delsq + rho * (z[1 - 1] * z[1 - 1] + z[2 - 1] * z[2 - 1]);
            c = rho * z[1 - 1] * z[1 - 1] * delsq;
            //
            //           B > ZERO, always
            //
            //           The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )
            //
            tau = two * c / (b + sqrt(abs(b * b - four * c)));
            //
            //           The following TAU is DSIGMA - D( 1 )
            //
            tau = tau / (d[1 - 1] + sqrt(d[1 - 1] * d[1 - 1] + tau));
            dsigma = d[1 - 1] + tau;
            delta[1 - 1] = -tau;
            delta[2 - 1] = del - tau;
            work[1 - 1] = two * d[1 - 1] + tau;
            work[2 - 1] = (d[1 - 1] + tau) + d[2 - 1];
            //           DELTA( 1 ) = -Z( 1 ) / TAU
            //           DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
        } else {
            b = -delsq + rho * (z[1 - 1] * z[1 - 1] + z[2 - 1] * z[2 - 1]);
            c = rho * z[2 - 1] * z[2 - 1] * delsq;
            //
            //           The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
            //
            if (b > zero) {
                tau = -two * c / (b + sqrt(b * b + four * c));
            } else {
                tau = (b - sqrt(b * b + four * c)) / two;
            }
            //
            //           The following TAU is DSIGMA - D( 2 )
            //
            tau = tau / (d[2 - 1] + sqrt(abs(d[2 - 1] * d[2 - 1] + tau)));
            dsigma = d[2 - 1] + tau;
            delta[1 - 1] = -(del + tau);
            delta[2 - 1] = -tau;
            work[1 - 1] = d[1 - 1] + tau + d[2 - 1];
            work[2 - 1] = two * d[2 - 1] + tau;
            //           DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            //           DELTA( 2 ) = -Z( 2 ) / TAU
        }
        //        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
        //        DELTA( 1 ) = DELTA( 1 ) / TEMP
        //        DELTA( 2 ) = DELTA( 2 ) / TEMP
    } else {
        //
        //        Now I=2
        //
        b = -delsq + rho * (z[1 - 1] * z[1 - 1] + z[2 - 1] * z[2 - 1]);
        c = rho * z[2 - 1] * z[2 - 1] * delsq;
        //
        //        The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
        //
        if (b > zero) {
            tau = (b + sqrt(b * b + four * c)) / two;
        } else {
            tau = two * c / (-b + sqrt(b * b + four * c));
        }
        //
        //        The following TAU is DSIGMA - D( 2 )
        //
        tau = tau / (d[2 - 1] + sqrt(d[2 - 1] * d[2 - 1] + tau));
        dsigma = d[2 - 1] + tau;
        delta[1 - 1] = -(del + tau);
        delta[2 - 1] = -tau;
        work[1 - 1] = d[1 - 1] + tau + d[2 - 1];
        work[2 - 1] = two * d[2 - 1] + tau;
        //        DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
        //        DELTA( 2 ) = -Z( 2 ) / TAU
        //        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
        //        DELTA( 1 ) = DELTA( 1 ) / TEMP
        //        DELTA( 2 ) = DELTA( 2 ) / TEMP
    }
    //
    //     End of Rlasd5
    //
}
