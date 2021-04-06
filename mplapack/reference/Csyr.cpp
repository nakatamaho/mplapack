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

void Csyr(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *a, INTEGER const lda) {
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
    //     .. Executable Statements ..
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
    } else if (lda < max((INTEGER)1, n)) {
        info = 7;
    }
    if (info != 0) {
        Mxerbla("Csyr  ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const COMPLEX zero = (0.0, 0.0);
    if ((n == 0) || (alpha == zero)) {
        return;
    }
    //
    //     Set the start poINTEGER in X if the increment is not unity.
    //
    INTEGER kx = 0;
    if (incx <= 0) {
        kx = 1 - (n - 1) * incx;
    } else if (incx != 1) {
        kx = 1;
    }
    //
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through the triangular part
    //     of A.
    //
    INTEGER j = 0;
    COMPLEX temp = 0.0;
    INTEGER i = 0;
    INTEGER jx = 0;
    INTEGER ix = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        Form  A  when A is stored in upper triangle.
        //
        if (incx == 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (x[j - 1] != zero) {
                    temp = alpha * x[j - 1];
                    for (i = 1; i <= j; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] += x[i - 1] * temp;
                    }
                }
            }
        } else {
            jx = kx;
            for (j = 1; j <= n; j = j + 1) {
                if (x[jx - 1] != zero) {
                    temp = alpha * x[jx - 1];
                    ix = kx;
                    for (i = 1; i <= j; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] += x[ix - 1] * temp;
                        ix += incx;
                    }
                }
                jx += incx;
            }
        }
    } else {
        //
        //        Form  A  when A is stored in lower triangle.
        //
        if (incx == 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (x[j - 1] != zero) {
                    temp = alpha * x[j - 1];
                    for (i = j; i <= n; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] += x[i - 1] * temp;
                    }
                }
            }
        } else {
            jx = kx;
            for (j = 1; j <= n; j = j + 1) {
                if (x[jx - 1] != zero) {
                    temp = alpha * x[jx - 1];
                    ix = jx;
                    for (i = j; i <= n; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] += x[ix - 1] * temp;
                        ix += incx;
                    }
                }
                jx += incx;
            }
        }
    }
    //
    //     End of Csyr
    //
}
