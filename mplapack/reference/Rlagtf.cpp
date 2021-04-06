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

void Rlagtf(INTEGER const n, REAL *a, REAL const lambda, REAL *b, REAL *c, REAL const tol, REAL *d, INTEGER *in, INTEGER &info) {
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
    // =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if (n < 0) {
        info = -1;
        Mxerbla("Rlagtf", -info);
        return;
    }
    //
    if (n == 0) {
        return;
    }
    //
    a[1 - 1] = a[1 - 1] - lambda;
    in[n - 1] = 0;
    const REAL zero = 0.0;
    if (n == 1) {
        if (a[1 - 1] == zero) {
            in[1 - 1] = 1;
        }
        return;
    }
    //
    REAL eps = Rlamch("Epsilon");
    //
    REAL tl = max(tol, eps);
    REAL scale1 = abs(a[1 - 1]) + abs(b[1 - 1]);
    INTEGER k = 0;
    REAL scale2 = 0.0;
    REAL piv1 = 0.0;
    REAL piv2 = 0.0;
    REAL mult = 0.0;
    REAL temp = 0.0;
    for (k = 1; k <= n - 1; k = k + 1) {
        a[(k + 1) - 1] = a[(k + 1) - 1] - lambda;
        scale2 = abs(c[k - 1]) + abs(a[(k + 1) - 1]);
        if (k < (n - 1)) {
            scale2 += abs(b[(k + 1) - 1]);
        }
        if (a[k - 1] == zero) {
            piv1 = zero;
        } else {
            piv1 = abs(a[k - 1]) / scale1;
        }
        if (c[k - 1] == zero) {
            in[k - 1] = 0;
            piv2 = zero;
            scale1 = scale2;
            if (k < (n - 1)) {
                d[k - 1] = zero;
            }
        } else {
            piv2 = abs(c[k - 1]) / scale2;
            if (piv2 <= piv1) {
                in[k - 1] = 0;
                scale1 = scale2;
                c[k - 1] = c[k - 1] / a[k - 1];
                a[(k + 1) - 1] = a[(k + 1) - 1] - c[k - 1] * b[k - 1];
                if (k < (n - 1)) {
                    d[k - 1] = zero;
                }
            } else {
                in[k - 1] = 1;
                mult = a[k - 1] / c[k - 1];
                a[k - 1] = c[k - 1];
                temp = a[(k + 1) - 1];
                a[(k + 1) - 1] = b[k - 1] - mult * temp;
                if (k < (n - 1)) {
                    d[k - 1] = b[(k + 1) - 1];
                    b[(k + 1) - 1] = -mult * d[k - 1];
                }
                b[k - 1] = temp;
                c[k - 1] = mult;
            }
        }
        if ((max(piv1, piv2) <= tl) && (in[n - 1] == 0)) {
            in[n - 1] = k;
        }
    }
    if ((abs(a[n - 1]) <= scale1 * tl) && (in[n - 1] == 0)) {
        in[n - 1] = n;
    }
    //
    //     End of Rlagtf
    //
}
