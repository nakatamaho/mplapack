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

void Clar2v(INTEGER const &n, COMPLEX *x, COMPLEX *y, COMPLEX *z, INTEGER const &incx, REAL *c, COMPLEX *s, INTEGER const &incc) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER ix = 1;
    INTEGER ic = 1;
    INTEGER i = 0;
    REAL xi = 0.0;
    REAL yi = 0.0;
    COMPLEX zi = 0.0;
    REAL zir = 0.0;
    REAL zii = 0.0;
    REAL ci = 0.0;
    COMPLEX si = 0.0;
    REAL sir = 0.0;
    REAL sii = 0.0;
    REAL t1r = 0.0;
    REAL t1i = 0.0;
    COMPLEX t2 = 0.0;
    COMPLEX t3 = 0.0;
    COMPLEX t4 = 0.0;
    REAL t5 = 0.0;
    REAL t6 = 0.0;
    for (i = 1; i <= n; i = i + 1) {
        xi = x[ix - 1].real();
        yi = y[ix - 1].real();
        zi = z[ix - 1];
        zir = zi.real();
        zii = zi.imag();
        ci = c[ic - 1];
        si = s[ic - 1];
        sir = si.real();
        sii = si.imag();
        t1r = sir * zir - sii * zii;
        t1i = sir * zii + sii * zir;
        t2 = ci * zi;
        t3 = t2 - conj(si) * xi;
        t4 = conj(t2) + si * yi;
        t5 = ci * xi + t1r;
        t6 = ci * yi - t1r;
        x[ix - 1] = ci * t5 + (sir * t4.real() + sii * t4.imag());
        y[ix - 1] = ci * t6 - (sir * t3.real() - sii * t3.imag());
        z[ix - 1] = ci * t3 + conj(si) * COMPLEX(t6, t1i);
        ix += incx;
        ic += incc;
    }
    //
    //     End of Clar2v
    //
}
