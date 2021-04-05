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

REAL Rasum(INTEGER const n, REAL *dx, INTEGER const incx) {
    REAL return_value = 0.0;
    //
    //  -- Reference BLAS level1 routine --
    //  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
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
    return_value = 0.0;
    REAL dtemp = 0.0;
    if (n <= 0 || incx <= 0) {
        return return_value;
    }
    INTEGER m = 0;
    INTEGER i = 0;
    INTEGER mp1 = 0;
    INTEGER nincx = 0;
    if (incx == 1) {
        //        code for increment equal to 1
        //
        //        clean-up loop
        //
        m = mod(n, 6);
        if (m != 0) {
            for (i = 1; i <= m; i = i + 1) {
                dtemp += abs(dx[i - 1]);
            }
            if (n < 6) {
                return_value = dtemp;
                return return_value;
            }
        }
        mp1 = m + 1;
        for (i = mp1; i <= n; i = i + 6) {
            dtemp += abs(dx[i - 1]) + abs(dx[(i + 1) - 1]) + abs(dx[(i + 2) - 1]) + abs(dx[(i + 3) - 1]) + abs(dx[(i + 4) - 1]) + abs(dx[(i + 5) - 1]);
        }
    } else {
        //
        //        code for increment not equal to 1
        //
        nincx = n * incx;
        for (i = 1; i <= nincx; i = i + incx) {
            dtemp += abs(dx[i - 1]);
        }
    }
    return_value = dtemp;
    return return_value;
}
