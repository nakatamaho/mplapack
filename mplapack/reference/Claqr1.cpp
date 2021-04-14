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

void Claqr1(INTEGER const n, COMPLEX *h, INTEGER const ldh, COMPLEX const s1, COMPLEX const s2, COMPLEX *v) {
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
    //  ================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX cdum = 0.0;
    abs1[cdum - 1] = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n != 2 && n != 3) {
        return;
    }
    //
    REAL s = 0.0;
    const REAL rzero = 0.0;
    const COMPLEX zero = (0.0, 0.0);
    COMPLEX h21s = 0.0;
    COMPLEX h31s = 0.0;
    if (n == 2) {
        s = abs1[(h[(1 - 1)] - s2) - 1] + abs1[h[(2 - 1)] - 1];
        if (s == rzero) {
            v[1 - 1] = zero;
            v[2 - 1] = zero;
        } else {
            h21s = h[(2 - 1)] / s;
            v[1 - 1] = h21s * h[(2 - 1) * ldh] + (h[(1 - 1)] - s1) * ((h[(1 - 1)] - s2) / s);
            v[2 - 1] = h21s * (h[(1 - 1)] + h[(2 - 1) + (2 - 1) * ldh] - s1 - s2);
        }
    } else {
        s = abs1[(h[(1 - 1)] - s2) - 1] + abs1[h[(2 - 1)] - 1] + abs1[h[(3 - 1)] - 1];
        if (s == zero) {
            v[1 - 1] = zero;
            v[2 - 1] = zero;
            v[3 - 1] = zero;
        } else {
            h21s = h[(2 - 1)] / s;
            h31s = h[(3 - 1)] / s;
            v[1 - 1] = (h[(1 - 1)] - s1) * ((h[(1 - 1)] - s2) / s) + h[(2 - 1) * ldh] * h21s + h[(3 - 1) * ldh] * h31s;
            v[2 - 1] = h21s * (h[(1 - 1)] + h[(2 - 1) + (2 - 1) * ldh] - s1 - s2) + h[(2 - 1) + (3 - 1) * ldh] * h31s;
            v[3 - 1] = h31s * (h[(1 - 1)] + h[(3 - 1) + (3 - 1) * ldh] - s1 - s2) + h21s * h[(3 - 1) + (2 - 1) * ldh];
        }
    }
}
