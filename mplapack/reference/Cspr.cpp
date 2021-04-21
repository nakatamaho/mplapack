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

void Cspr(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *ap) {
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
    }
    if (info != 0) {
        Mxerbla("Cspr  ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    if ((n == 0) || (alpha == zero)) {
        return;
    }
    //
    //     Set the start point in X if the increment is not unity.
    //
    INTEGER kx = 0;
    if (incx <= 0) {
        kx = 1 - (n - 1) * incx;
    } else if (incx != 1) {
        kx = 1;
    }
    //
    //     Start the operations. In this version the elements of the array AP
    //     are accessed sequentially with one pass through AP.
    //
    INTEGER kk = 1;
    INTEGER j = 0;
    COMPLEX temp = 0.0;
    INTEGER k = 0;
    INTEGER i = 0;
    INTEGER jx = 0;
    INTEGER ix = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        Form  A  when upper triangle is stored in AP.
        //
        if (incx == 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (x[j - 1] != zero) {
                    temp = alpha * x[j - 1];
                    k = kk;
                    for (i = 1; i <= j - 1; i = i + 1) {
                        ap[k - 1] += x[i - 1] * temp;
                        k++;
                    }
                    ap[(kk + j - 1) - 1] += x[j - 1] * temp;
                } else {
                    ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1];
                }
                kk += j;
            }
        } else {
            jx = kx;
            for (j = 1; j <= n; j = j + 1) {
                if (x[jx - 1] != zero) {
                    temp = alpha * x[jx - 1];
                    ix = kx;
                    for (k = kk; k <= kk + j - 2; k = k + 1) {
                        ap[k - 1] += x[ix - 1] * temp;
                        ix += incx;
                    }
                    ap[(kk + j - 1) - 1] += x[jx - 1] * temp;
                } else {
                    ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1];
                }
                jx += incx;
                kk += j;
            }
        }
    } else {
        //
        //        Form  A  when lower triangle is stored in AP.
        //
        if (incx == 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (x[j - 1] != zero) {
                    temp = alpha * x[j - 1];
                    ap[kk - 1] += temp * x[j - 1];
                    k = kk + 1;
                    for (i = j + 1; i <= n; i = i + 1) {
                        ap[k - 1] += x[i - 1] * temp;
                        k++;
                    }
                } else {
                    ap[kk - 1] = ap[kk - 1];
                }
                kk += n - j + 1;
            }
        } else {
            jx = kx;
            for (j = 1; j <= n; j = j + 1) {
                if (x[jx - 1] != zero) {
                    temp = alpha * x[jx - 1];
                    ap[kk - 1] += temp * x[jx - 1];
                    ix = jx;
                    for (k = kk + 1; k <= kk + n - j; k = k + 1) {
                        ix += incx;
                        ap[k - 1] += x[ix - 1] * temp;
                    }
                } else {
                    ap[kk - 1] = ap[kk - 1];
                }
                jx += incx;
                kk += n - j + 1;
            }
        }
    }
    //
    //     End of Cspr
    //
}
