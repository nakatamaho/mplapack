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

void Rlar2v(INTEGER const n, REAL *x, REAL *y, REAL *z, INTEGER const incx, REAL *c, REAL *s, INTEGER const incc) {
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
    //     .. Local Scalars ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER ix = 1;
    INTEGER ic = 1;
    INTEGER i = 0;
    REAL xi = 0.0;
    REAL yi = 0.0;
    REAL zi = 0.0;
    REAL ci = 0.0;
    REAL si = 0.0;
    REAL t1 = 0.0;
    REAL t2 = 0.0;
    REAL t3 = 0.0;
    REAL t4 = 0.0;
    REAL t5 = 0.0;
    REAL t6 = 0.0;
    for (i = 1; i <= n; i = i + 1) {
        xi = x[ix - 1];
        yi = y[ix - 1];
        zi = z[ix - 1];
        ci = c[ic - 1];
        si = s[ic - 1];
        t1 = si * zi;
        t2 = ci * zi;
        t3 = t2 - si * xi;
        t4 = t2 + si * yi;
        t5 = ci * xi + t1;
        t6 = ci * yi - t1;
        x[ix - 1] = ci * t5 + si * t4;
        y[ix - 1] = ci * t6 - si * t3;
        z[ix - 1] = ci * t4 - si * t5;
        ix += incx;
        ic += incc;
    }
    //
    //     End of Rlar2v
    //
}
