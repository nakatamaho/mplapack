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

void Rlaqr1(INTEGER const &n, REAL *h, INTEGER const &ldh, REAL const &sr1, REAL const &si1, REAL const &sr2, REAL const &si2, REAL *v) {
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
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n != 2 && n != 3) {
        return;
    }
    //
    REAL s = 0.0;
    const REAL zero = 0.0;
    REAL h21s = 0.0;
    REAL h31s = 0.0;
    if (n == 2) {
        s = abs(h[(1 - 1)] - sr2) + abs(si2) + abs(h[(2 - 1)]);
        if (s == zero) {
            v[1 - 1] = zero;
            v[2 - 1] = zero;
        } else {
            h21s = h[(2 - 1)] / s;
            v[1 - 1] = h21s * h[(2 - 1) * ldh] + (h[(1 - 1)] - sr1) * ((h[(1 - 1)] - sr2) / s) - si1 * (si2 / s);
            v[2 - 1] = h21s * (h[(1 - 1)] + h[(2 - 1) + (2 - 1) * ldh] - sr1 - sr2);
        }
    } else {
        s = abs(h[(1 - 1)] - sr2) + abs(si2) + abs(h[(2 - 1)]) + abs(h[(3 - 1)]);
        if (s == zero) {
            v[1 - 1] = zero;
            v[2 - 1] = zero;
            v[3 - 1] = zero;
        } else {
            h21s = h[(2 - 1)] / s;
            h31s = h[(3 - 1)] / s;
            v[1 - 1] = (h[(1 - 1)] - sr1) * ((h[(1 - 1)] - sr2) / s) - si1 * (si2 / s) + h[(2 - 1) * ldh] * h21s + h[(3 - 1) * ldh] * h31s;
            v[2 - 1] = h21s * (h[(1 - 1)] + h[(2 - 1) + (2 - 1) * ldh] - sr1 - sr2) + h[(2 - 1) + (3 - 1) * ldh] * h31s;
            v[3 - 1] = h31s * (h[(1 - 1)] + h[(3 - 1) + (3 - 1) * ldh] - sr1 - sr2) + h21s * h[(3 - 1) + (2 - 1) * ldh];
        }
    }
}
