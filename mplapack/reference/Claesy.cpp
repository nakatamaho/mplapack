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

void Claesy(COMPLEX const a, COMPLEX const b, COMPLEX const c, COMPLEX &rt1, COMPLEX &rt2, COMPLEX &evscal, COMPLEX &cs1, COMPLEX &sn1) {
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
    //     Special case:  The matrix is actually diagonal.
    //     To avoid divide by zero later, we treat this case separately.
    //
    const REAL zero = 0.0;
    COMPLEX tmp = 0.0;
    const REAL one = 1.0;
    const REAL half = 0.5e0;
    COMPLEX s = 0.0;
    COMPLEX t = 0.0;
    REAL babs = 0.0;
    REAL tabs = 0.0;
    REAL z = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    REAL evnorm = 0.0;
    const REAL thresh = 0.1e0;
    if (abs(b) == zero) {
        rt1 = a;
        rt2 = c;
        if (abs(rt1) < abs(rt2)) {
            tmp = rt1;
            rt1 = rt2;
            rt2 = tmp;
            cs1 = zero;
            sn1 = one;
        } else {
            cs1 = one;
            sn1 = zero;
        }
    } else {
        //
        //        Compute the eigenvalues and eigenvectors.
        //        The characteristic equation is
        //           lambda **2 - (A+C) lambda + (A*C - B*B)
        //        and we solve it using the quadratic formula.
        //
        s = (a + c) * half;
        t = (a - c) * half;
        //
        //        Take the square root carefully to avoid over/under flow.
        //
        babs = abs(b);
        tabs = abs(t);
        z = max(babs, tabs);
        if (z > zero) {
            t = z * sqrt((t / z) * (t / z) + (b / z) * (b / z));
        }
        //
        //        Compute the two eigenvalues.  RT1 and RT2 are exchanged
        //        if necessary so that RT1 will have the greater magnitude.
        //
        rt1 = s + t;
        rt2 = s - t;
        if (abs(rt1) < abs(rt2)) {
            tmp = rt1;
            rt1 = rt2;
            rt2 = tmp;
        }
        //
        //        Choose CS1 = 1 and SN1 to satisfy the first equation, then
        //        scale the components of this eigenvector so that the matrix
        //        of eigenvectors X satisfies  X * X**T = I .  (No scaling is
        //        done if the norm of the eigenvalue matrix is less than THRESH.)
        //
        sn1 = (rt1 - a) / b;
        tabs = abs(sn1);
        if (tabs > one) {
            t = tabs * sqrt((one / tabs) * (one / tabs) + (sn1 / tabs) * (sn1 / tabs));
        } else {
            t = sqrt(cone + sn1 * sn1);
        }
        evnorm = abs(t);
        if (evnorm >= thresh) {
            evscal = cone / t;
            cs1 = evscal;
            sn1 = sn1 * evscal;
        } else {
            evscal = zero;
        }
    }
    //
    //     End of Claesy
    //
}
