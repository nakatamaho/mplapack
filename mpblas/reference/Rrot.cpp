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

void Rrot(INTEGER const &n, REAL *dx, INTEGER const &incx, REAL *dy, INTEGER const &incy, REAL const &c, REAL const &s) {
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
    if (n <= 0) {
        return;
    }
    INTEGER i = 0;
    REAL dtemp = 0.0;
    INTEGER ix = 0;
    INTEGER iy = 0;
    if (incx == 1 && incy == 1) {
        //
        //       code for both increments equal to 1
        //
        for (i = 1; i <= n; i = i + 1) {
            dtemp = c * dx[i - 1] + s * dy[i - 1];
            dy[i - 1] = c * dy[i - 1] - s * dx[i - 1];
            dx[i - 1] = dtemp;
        }
    } else {
        //
        //       code for unequal increments or equal increments not equal
        //         to 1
        //
        ix = 1;
        iy = 1;
        if (incx < 0) {
            ix = (-n + 1) * incx + 1;
        }
        if (incy < 0) {
            iy = (-n + 1) * incy + 1;
        }
        for (i = 1; i <= n; i = i + 1) {
            dtemp = c * dx[ix - 1] + s * dy[iy - 1];
            dy[iy - 1] = c * dy[iy - 1] - s * dx[ix - 1];
            dx[ix - 1] = dtemp;
            ix += incx;
            iy += incy;
        }
    }
}
