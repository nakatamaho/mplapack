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

REAL Cla_gerpvgrw(INTEGER const n, INTEGER const ncols, COMPLEX *a, INTEGER const lda, COMPLEX *af, INTEGER const ldaf) {
    REAL return_value = 0.0;
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
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function Definitions ..
    COMPLEX zdum = 0.0;
    abs1(zdum) = abs(zdum.real()) + abs(zdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    REAL rpvgrw = 1.0;
    //
    INTEGER j = 0;
    REAL amax = 0.0;
    REAL umax = 0.0;
    INTEGER i = 0;
    for (j = 1; j <= ncols; j = j + 1) {
        amax = 0.0;
        umax = 0.0;
        for (i = 1; i <= n; i = i + 1) {
            amax = max(abs1(a[(i - 1) + (j - 1) * lda]), amax);
        }
        for (i = 1; i <= j; i = i + 1) {
            umax = max(abs1(af[(i - 1) + (j - 1) * ldaf]), umax);
        }
        if (umax != 0.0) {
            rpvgrw = min(amax / umax, rpvgrw);
        }
    }
    return_value = rpvgrw;
    return return_value;
}
