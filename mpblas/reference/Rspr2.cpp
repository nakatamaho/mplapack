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

void Rspr2(const char *uplo, INTEGER const &n, REAL const &alpha, REAL *x, INTEGER const &incx, REAL *y, INTEGER const &incy, REAL *ap) {
    //
    //  -- Reference BLAS level2 routine --
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
        info = 1;
    } else if (n < 0) {
        info = 2;
    } else if (incx == 0) {
        info = 5;
    } else if (incy == 0) {
        info = 7;
    }
    if (info != 0) {
        Mxerbla("Rspr2 ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const REAL zero = 0.0;
    if ((n == 0) || (alpha == zero)) {
        return;
    }
    //
    //     Set up the start poINTEGERs in X and Y if the increments are not both
    //     unity.
    //
    INTEGER kx = 0;
    INTEGER ky = 0;
    INTEGER jx = 0;
    INTEGER jy = 0;
    if ((incx != 1) || (incy != 1)) {
        if (incx > 0) {
            kx = 1;
        } else {
            kx = 1 - (n - 1) * incx;
        }
        if (incy > 0) {
            ky = 1;
        } else {
            ky = 1 - (n - 1) * incy;
        }
        jx = kx;
        jy = ky;
    }
    //
    //     Start the operations. In this version the elements of the array AP
    //     are accessed sequentially with one pass through AP.
    //
    INTEGER kk = 1;
    INTEGER j = 0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    INTEGER k = 0;
    INTEGER i = 0;
    INTEGER ix = 0;
    INTEGER iy = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        Form  A  when upper triangle is stored in AP.
        //
        if ((incx == 1) && (incy == 1)) {
            for (j = 1; j <= n; j = j + 1) {
                if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                    temp1 = alpha * y[j - 1];
                    temp2 = alpha * x[j - 1];
                    k = kk;
                    for (i = 1; i <= j; i = i + 1) {
                        ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                        k++;
                    }
                }
                kk += j;
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                if ((x[jx - 1] != zero) || (y[jy - 1] != zero)) {
                    temp1 = alpha * y[jy - 1];
                    temp2 = alpha * x[jx - 1];
                    ix = kx;
                    iy = ky;
                    for (k = kk; k <= kk + j - 1; k = k + 1) {
                        ap[k - 1] += x[ix - 1] * temp1 + y[iy - 1] * temp2;
                        ix += incx;
                        iy += incy;
                    }
                }
                jx += incx;
                jy += incy;
                kk += j;
            }
        }
    } else {
        //
        //        Form  A  when lower triangle is stored in AP.
        //
        if ((incx == 1) && (incy == 1)) {
            for (j = 1; j <= n; j = j + 1) {
                if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                    temp1 = alpha * y[j - 1];
                    temp2 = alpha * x[j - 1];
                    k = kk;
                    for (i = j; i <= n; i = i + 1) {
                        ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                        k++;
                    }
                }
                kk += n - j + 1;
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                if ((x[jx - 1] != zero) || (y[jy - 1] != zero)) {
                    temp1 = alpha * y[jy - 1];
                    temp2 = alpha * x[jx - 1];
                    ix = jx;
                    iy = jy;
                    for (k = kk; k <= kk + n - j; k = k + 1) {
                        ap[k - 1] += x[ix - 1] * temp1 + y[iy - 1] * temp2;
                        ix += incx;
                        iy += incy;
                    }
                }
                jx += incx;
                jy += incy;
                kk += n - j + 1;
            }
        }
    }
    //
    //     End of Rspr2 .
    //
}
