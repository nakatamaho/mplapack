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

void Rlargv(INTEGER const &n, REAL *x, INTEGER const &incx, REAL *y, INTEGER const &incy, REAL *c, INTEGER const &incc) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER ix = 1;
    INTEGER iy = 1;
    INTEGER ic = 1;
    INTEGER i = 0;
    REAL f = 0.0;
    REAL g = 0.0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL t = 0.0;
    REAL tt = 0.0;
    for (i = 1; i <= n; i = i + 1) {
        f = x[ix - 1];
        g = y[iy - 1];
        if (g == zero) {
            c[ic - 1] = one;
        } else if (f == zero) {
            c[ic - 1] = zero;
            y[iy - 1] = one;
            x[ix - 1] = g;
        } else if (abs(f) > abs(g)) {
            t = g / f;
            tt = sqrt(one + t * t);
            c[ic - 1] = one / tt;
            y[iy - 1] = t * c[ic - 1];
            x[ix - 1] = f * tt;
        } else {
            t = f / g;
            tt = sqrt(one + t * t);
            y[iy - 1] = one / tt;
            c[ic - 1] = t * y[iy - 1];
            x[ix - 1] = g * tt;
        }
        ic += incc;
        iy += incy;
        ix += incx;
    }
    //
    //     End of Rlargv
    //
}
