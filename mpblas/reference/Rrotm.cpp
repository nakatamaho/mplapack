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

void Rrotm(INTEGER const n, REAL *dx, INTEGER const incx, REAL *dy, INTEGER const incy, REAL *dparam) {
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
    REAL zero = 0.0;
    REAL two = 2.0;
    //     ..
    //
    REAL dflag = dparam[1 - 1];
    if (n <= 0 || (dflag + two == zero)) {
        return;
    }
    INTEGER nsteps = 0;
    REAL dh11 = 0.0;
    REAL dh12 = 0.0;
    REAL dh21 = 0.0;
    REAL dh22 = 0.0;
    INTEGER i = 0;
    REAL w = 0.0;
    REAL z = 0.0;
    INTEGER kx = 0;
    INTEGER ky = 0;
    if (incx == incy && incx > 0) {
        //
        nsteps = n * incx;
        if (dflag < zero) {
            dh11 = dparam[2 - 1];
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= nsteps; i = i + incx) {
                w = dx[i - 1];
                z = dy[i - 1];
                dx[i - 1] = w * dh11 + z * dh12;
                dy[i - 1] = w * dh21 + z * dh22;
            }
        } else if (dflag == zero) {
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            for (i = 1; i <= nsteps; i = i + incx) {
                w = dx[i - 1];
                z = dy[i - 1];
                dx[i - 1] = w + z * dh12;
                dy[i - 1] = w * dh21 + z;
            }
        } else {
            dh11 = dparam[2 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= nsteps; i = i + incx) {
                w = dx[i - 1];
                z = dy[i - 1];
                dx[i - 1] = w * dh11 + z;
                dy[i - 1] = -w + dh22 * z;
            }
        }
    } else {
        kx = 1;
        ky = 1;
        if (incx < 0) {
            kx = 1 + (1 - n) * incx;
        }
        if (incy < 0) {
            ky = 1 + (1 - n) * incy;
        }
        //
        if (dflag < zero) {
            dh11 = dparam[2 - 1];
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= n; i = i + 1) {
                w = dx[kx - 1];
                z = dy[ky - 1];
                dx[kx - 1] = w * dh11 + z * dh12;
                dy[ky - 1] = w * dh21 + z * dh22;
                kx += incx;
                ky += incy;
            }
        } else if (dflag == zero) {
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            for (i = 1; i <= n; i = i + 1) {
                w = dx[kx - 1];
                z = dy[ky - 1];
                dx[kx - 1] = w + z * dh12;
                dy[ky - 1] = w * dh21 + z;
                kx += incx;
                ky += incy;
            }
        } else {
            dh11 = dparam[2 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= n; i = i + 1) {
                w = dx[kx - 1];
                z = dy[ky - 1];
                dx[kx - 1] = w * dh11 + z;
                dy[ky - 1] = -w + dh22 * z;
                kx += incx;
                ky += incy;
            }
        }
    }
}
