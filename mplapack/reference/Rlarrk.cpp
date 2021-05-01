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

void Rlarrk(INTEGER const n, INTEGER const iw, REAL const gl, REAL const gu, REAL *d, REAL *e2, REAL const pivmin, REAL const reltol, REAL &w, REAL &werr, INTEGER &info) {
    REAL eps = 0.0;
    REAL tnorm = 0.0;
    REAL rtoli = 0.0;
    const REAL two = 2.0;
    const REAL fudge = two;
    REAL atoli = 0.0;
    INTEGER itmax = 0;
    REAL left = 0.0;
    REAL right = 0.0;
    INTEGER it = 0;
    REAL tmp1 = 0.0;
    REAL tmp2 = 0.0;
    const REAL half = 0.5e0;
    REAL mid = 0.0;
    INTEGER negcnt = 0;
    const REAL zero = 0.0;
    INTEGER i = 0;
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        info = 0;
        return;
    }
    //
    //     Get machine constants
    eps = Rlamch("P");
    //
    tnorm = max(abs(gl), abs(gu));
    rtoli = reltol;
    atoli = fudge * two * pivmin;
    //
    itmax = castREAL((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
    //
    info = -1;
    //
    left = gl - fudge * tnorm * eps * n - fudge * two * pivmin;
    right = gu + fudge * tnorm * eps * n + fudge * two * pivmin;
    it = 0;
//
statement_10:
    //
    //     Check if interval converged or maximum number of iterations reached
    //
    tmp1 = abs(right - left);
    tmp2 = max(abs(right), abs(left));
    if (tmp1 < max({atoli, pivmin, rtoli * tmp2})) {
        info = 0;
        goto statement_30;
    }
    if (it > itmax) {
        goto statement_30;
    }
    //
    //     Count number of negative pivots for mid-point
    //
    it++;
    mid = half * (left + right);
    negcnt = 0;
    tmp1 = d[1 - 1] - mid;
    if (abs(tmp1) < pivmin) {
        tmp1 = -pivmin;
    }
    if (tmp1 <= zero) {
        negcnt++;
    }
    //
    for (i = 2; i <= n; i = i + 1) {
        tmp1 = d[i - 1] - e2[(i - 1) - 1] / tmp1 - mid;
        if (abs(tmp1) < pivmin) {
            tmp1 = -pivmin;
        }
        if (tmp1 <= zero) {
            negcnt++;
        }
    }
    //
    if (negcnt >= iw) {
        right = mid;
    } else {
        left = mid;
    }
    goto statement_10;
//
statement_30:
    //
    //     Converged or maximum number of iterations reached
    //
    w = half * (left + right);
    werr = half * abs(right - left);
    //
    //     End of Rlarrk
    //
}
