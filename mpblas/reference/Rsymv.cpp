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

void Rsymv(const char *uplo, INTEGER const &n, REAL const &alpha, REAL *a, INTEGER const &lda, REAL *x, INTEGER const &incx, REAL const &beta, REAL *y, INTEGER const &incy) {
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
    } else if (lda < max((INTEGER)1, n)) {
        info = 5;
    } else if (incx == 0) {
        info = 7;
    } else if (incy == 0) {
        info = 10;
    }
    if (info != 0) {
        Mxerbla("Rsymv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if ((n == 0) || ((alpha == zero) && (beta == one))) {
        return;
    }
    //
    //     Set up the start poINTEGERs in  X  and  Y.
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
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through the triangular part
    //     of A.
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
    INTEGER j = 0;
    REAL temp1 = 0.0;
    REAL temp2 = 0.0;
    INTEGER jx = 0;
    INTEGER jy = 0;
    INTEGER ix = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        Form  y  when A is stored in upper triangle.
        //
        if ((incx == 1) && (incy == 1)) {
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[j - 1];
                temp2 = zero;
                for (i = 1; i <= j - 1; i = i + 1) {
                    y[i - 1] += temp1 * a[(i - 1) + (j - 1) * lda];
                    temp2 += a[(i - 1) + (j - 1) * lda] * x[i - 1];
                }
                y[j - 1] += temp1 * a[(j - 1) + (j - 1) * lda] + alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[jx - 1];
                temp2 = zero;
                ix = kx;
                iy = ky;
                for (i = 1; i <= j - 1; i = i + 1) {
                    y[iy - 1] += temp1 * a[(i - 1) + (j - 1) * lda];
                    temp2 += a[(i - 1) + (j - 1) * lda] * x[ix - 1];
                    ix += incx;
                    iy += incy;
                }
                y[jy - 1] += temp1 * a[(j - 1) + (j - 1) * lda] + alpha * temp2;
                jx += incx;
                jy += incy;
            }
        }
    } else {
        //
        //        Form  y  when A is stored in lower triangle.
        //
        if ((incx == 1) && (incy == 1)) {
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[j - 1];
                temp2 = zero;
                y[j - 1] += temp1 * a[(j - 1) + (j - 1) * lda];
                for (i = j + 1; i <= n; i = i + 1) {
                    y[i - 1] += temp1 * a[(i - 1) + (j - 1) * lda];
                    temp2 += a[(i - 1) + (j - 1) * lda] * x[i - 1];
                }
                y[j - 1] += alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            for (j = 1; j <= n; j = j + 1) {
                temp1 = alpha * x[jx - 1];
                temp2 = zero;
                y[jy - 1] += temp1 * a[(j - 1) + (j - 1) * lda];
                ix = jx;
                iy = jy;
                for (i = j + 1; i <= n; i = i + 1) {
                    ix += incx;
                    iy += incy;
                    y[iy - 1] += temp1 * a[(i - 1) + (j - 1) * lda];
                    temp2 += a[(i - 1) + (j - 1) * lda] * x[ix - 1];
                }
                y[jy - 1] += alpha * temp2;
                jx += incx;
                jy += incy;
            }
        }
    }
    //
    //     End of Rsymv .
    //
}
