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

void Rlaswp(INTEGER const n, REAL *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx) {
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
    //     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
    //     K1 through K2.
    //
    INTEGER ix0 = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER inc = 0;
    if (incx > 0) {
        ix0 = k1;
        i1 = k1;
        i2 = k2;
        inc = 1;
    } else if (incx < 0) {
        ix0 = k1 + (k1 - k2) * incx;
        i1 = k2;
        i2 = k1;
        inc = -1;
    } else {
        return;
    }
    //
    INTEGER n32 = (n / 32) * 32;
    INTEGER j = 0;
    INTEGER ix = 0;
    INTEGER i = 0;
    INTEGER ip = 0;
    INTEGER k = 0;
    REAL temp = 0.0;
    if (n32 != 0) {
        for (j = 1; j <= n32; j = j + 32) {
            ix = ix0;
            for (i = i1; inc > 0 ? i <= i2 : i >= i2; i = i + inc) {
                ip = ipiv[ix - 1];
                if (ip != i) {
                    for (k = j; k <= j + 31; k = k + 1) {
                        temp = a[(i - 1) + (k - 1) * lda];
                        a[(i - 1) + (k - 1) * lda] = a[(ip - 1) + (k - 1) * lda];
                        a[(ip - 1) + (k - 1) * lda] = temp;
                    }
                }
                ix += incx;
            }
        }
    }
    if (n32 != n) {
        n32++;
        ix = ix0;
        for (i = i1; inc > 0 ? i <= i2 : i >= i2; i = i + inc) {
            ip = ipiv[ix - 1];
            if (ip != i) {
                for (k = n32; k <= n; k = k + 1) {
                    temp = a[(i - 1) + (k - 1) * lda];
                    a[(i - 1) + (k - 1) * lda] = a[(ip - 1) + (k - 1) * lda];
                    a[(ip - 1) + (k - 1) * lda] = temp;
                }
            }
            ix += incx;
        }
    }
    //
    //     End of Rlaswp
    //
}
