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

void Clapll(INTEGER const &n, COMPLEX *x, INTEGER const &incx, COMPLEX *y, INTEGER const &incy, REAL &ssmin) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (n <= 1) {
        ssmin = zero;
        return;
    }
    //
    //     Compute the QR factorization of the N-by-2 matrix ( X Y )
    //
    COMPLEX tau = 0.0;
    Clarfg(n, x[1 - 1], x[(1 + incx) - 1], incx, tau);
    COMPLEX a11 = x[1 - 1];
    const COMPLEX cone = (1.0, 0.0);
    x[1 - 1] = cone;
    //
    COMPLEX c = -conj(tau) * Cdotc[(n - 1) + (x - 1) * ldCdotc];
    Caxpy(n, c, x, incx, y, incy);
    //
    Clarfg(n - 1, y[(1 + incy) - 1], y[(1 + 2 * incy) - 1], incy, tau);
    //
    COMPLEX a12 = y[1 - 1];
    COMPLEX a22 = y[(1 + incy) - 1];
    //
    //     Compute the SVD of 2-by-2 Upper triangular matrix.
    //
    REAL ssmax = 0.0;
    Rlas2(abs(a11), abs(a12), abs(a22), ssmin, ssmax);
    //
    //     End of Clapll
    //
}
