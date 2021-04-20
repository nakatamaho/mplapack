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

void Chbmv(const char *uplo, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
        info = 1;
    } else if (n < 0) {
        info = 2;
    } else if (k < 0) {
        info = 3;
    } else if (lda < (k + 1)) {
        info = 6;
    } else if (incx == 0) {
        info = 8;
    } else if (incy == 0) {
        info = 11;
    }
    if (info != 0) {
        Mxerbla("Chbmv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    if ((n == 0) || ((alpha == zero) && (beta == one))) {
        return;
    }
    //
    //     Set up the start points in  X  and  Y.
    //
    INTEGER kx = 0;
    if (incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (n - 1) * incx;
    }
    INTEGER ky = 0;
    if (incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (n - 1) * incy;
    }
    //
    //     Start the operations. In this version the elements of the array A
    //     are accessed sequentially with one pass through A.
    //
    //     First form  y := beta*y.
    //
    INTEGER i = 0;
    INTEGER iy = 0;
    if (beta != one) {
        if (incy == 1) {
            if (beta == zero) {
                for (i = 1; i <= n; i = i + 1) {
                    y[i - 1] = zero;
                }
            } else {
                for (i = 1; i <= n; i = i + 1) {
                    y[i - 1] = beta * y[i - 1];
                }
            }
        } else {
            iy = ky;
            if (beta == zero) {
                for (i = 1; i <= n; i = i + 1) {
                    y[iy - 1] = zero;
                    iy += incy;
                }
            } else {
                for (i = 1; i <= n; i = i + 1) {
                    y[iy - 1] = beta * y[iy - 1];
                    iy += incy;
                }
            }
        }
    }
    if (alpha == zero) {
        return;
    }
    INTEGER kplus1 = 0;
    INTEGER j = 0;
    COMPLEX temp1 = 0.0;
    COMPLEX temp2 = 0.0;
    INTEGER l = 0;
    INTEGER jx = 0;
    INTEGER jy = 0;
    INTEGER ix = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        Form  y  when upper triangle of A is stored.
        //
        kplus1 = k + 1;
        if ((incx == 1) && (incy == 1)) {
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[j - 1];
                temp2 = zero;
                l = kplus1 - j;
                for (i = max((INTEGER)1, j - k); i <= j - 1; i = i + 1) {
                    y[i - 1] += temp1 * a[((l + i) - 1) + (j - 1) * lda];
                    temp2 += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[i - 1];
                }
                y[j - 1] += temp1 * a[(kplus1 - 1) + (j - 1) * lda].real() + alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[jx - 1];
                temp2 = zero;
                ix = kx;
                iy = ky;
                l = kplus1 - j;
                for (i = max((INTEGER)1, j - k); i <= j - 1; i = i + 1) {
                    y[iy - 1] += temp1 * a[((l + i) - 1) + (j - 1) * lda];
                    temp2 += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[ix - 1];
                    ix += incx;
                    iy += incy;
                }
                y[jy - 1] += temp1 * a[(kplus1 - 1) + (j - 1) * lda].real() + alpha * temp2;
                jx += incx;
                jy += incy;
                if (j > k) {
                    kx += incx;
                    ky += incy;
                }
            }
        }
    } else {
        //
        //        Form  y  when lower triangle of A is stored.
        //
        if ((incx == 1) && (incy == 1)) {
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[j - 1];
                temp2 = zero;
                y[j - 1] += temp1 * a[(j - 1) * lda].real();
                l = 1 - j;
                for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                    y[i - 1] += temp1 * a[((l + i) - 1) + (j - 1) * lda];
                    temp2 += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[i - 1];
                }
                y[j - 1] += alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[jx - 1];
                temp2 = zero;
                y[jy - 1] += temp1 * a[(j - 1) * lda].real();
                l = 1 - j;
                ix = jx;
                iy = jy;
                for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                    ix += incx;
                    iy += incy;
                    y[iy - 1] += temp1 * a[((l + i) - 1) + (j - 1) * lda];
                    temp2 += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[ix - 1];
                }
                y[jy - 1] += alpha * temp2;
                jx += incx;
                jy += incy;
            }
        }
    }
    //
    //     End of Chbmv .
    //
}
