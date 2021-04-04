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

void Clacrt(INTEGER const &n, COMPLEX *cx, INTEGER const &incx, COMPLEX *cy, INTEGER const &incy, COMPLEX const &c, COMPLEX const &s) {
    INTEGER ix = 0;
    INTEGER iy = 0;
    INTEGER i = 0;
    COMPLEX ctemp = 0.0;
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
    // =====================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. Executable Statements ..
    //
    if (n <= 0) {
        return;
    }
    if (incx == 1 && incy == 1) {
        goto statement_20;
    }
    //
    //     Code for unequal increments or equal increments not equal to 1
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
        ctemp = c * cx[ix - 1] + s * cy[iy - 1];
        cy[iy - 1] = c * cy[iy - 1] - s * cx[ix - 1];
        cx[ix - 1] = ctemp;
        ix += incx;
        iy += incy;
    }
    return;
//
//     Code for both increments equal to 1
//
statement_20:
    for (i = 1; i <= n; i = i + 1) {
        ctemp = c * cx[i - 1] + s * cy[i - 1];
        cy[i - 1] = c * cy[i - 1] - s * cx[i - 1];
        cx[i - 1] = ctemp;
    }
}
