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

void Cheswapr(const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER const &lda, INTEGER const &i1, INTEGER const &i2) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //
    //  =====================================================================
    //
    //     ..
    //     .. Local Scalars ..
    //
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    bool upper = Mlsame(uplo, "U");
    COMPLEX tmp = 0.0;
    INTEGER i = 0;
    if (upper) {
        //
        //         UPPER
        //         first swap
        //          - swap column I1 and I2 from I1 to I1-1
        Cswap(i1 - 1, a[(i1 - 1) * lda], 1, a[(i2 - 1) * lda], 1);
        //
        //          second swap :
        //          - swap A(I1,I1) and A(I2,I2)
        //          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
        //          - swap A(I2,I1) and A(I1,I2)
        //
        tmp = a[(i1 - 1) + (i1 - 1) * lda];
        a[(i1 - 1) + (i1 - 1) * lda] = a[(i2 - 1) + (i2 - 1) * lda];
        a[(i2 - 1) + (i2 - 1) * lda] = tmp;
        //
        for (i = 1; i <= i2 - i1 - 1; i = i + 1) {
            tmp = a[(i1 - 1) + ((i1 + i) - 1) * lda];
            a[(i1 - 1) + ((i1 + i) - 1) * lda] = conj(a[((i1 + i) - 1) + (i2 - 1) * lda]);
            a[((i1 + i) - 1) + (i2 - 1) * lda] = conj(tmp);
        }
        //
        a[(i1 - 1) + (i2 - 1) * lda] = conj(a[(i1 - 1) + (i2 - 1) * lda]);
        //
        //          third swap
        //          - swap row I1 and I2 from I2+1 to N
        for (i = i2 + 1; i <= n; i = i + 1) {
            tmp = a[(i1 - 1) + (i - 1) * lda];
            a[(i1 - 1) + (i - 1) * lda] = a[(i2 - 1) + (i - 1) * lda];
            a[(i2 - 1) + (i - 1) * lda] = tmp;
        }
        //
    } else {
        //
        //         LOWER
        //         first swap
        //          - swap row I1 and I2 from 1 to I1-1
        Cswap(i1 - 1, a[(i1 - 1)], lda, a[(i2 - 1)], lda);
        //
        //         second swap :
        //          - swap A(I1,I1) and A(I2,I2)
        //          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
        //          - swap A(I2,I1) and A(I1,I2)
        //
        tmp = a[(i1 - 1) + (i1 - 1) * lda];
        a[(i1 - 1) + (i1 - 1) * lda] = a[(i2 - 1) + (i2 - 1) * lda];
        a[(i2 - 1) + (i2 - 1) * lda] = tmp;
        //
        for (i = 1; i <= i2 - i1 - 1; i = i + 1) {
            tmp = a[((i1 + i) - 1) + (i1 - 1) * lda];
            a[((i1 + i) - 1) + (i1 - 1) * lda] = conj(a[(i2 - 1) + ((i1 + i) - 1) * lda]);
            a[(i2 - 1) + ((i1 + i) - 1) * lda] = conj(tmp);
        }
        //
        a[(i2 - 1) + (i1 - 1) * lda] = conj(a[(i2 - 1) + (i1 - 1) * lda]);
        //
        //         third swap
        //          - swap col I1 and I2 from I2+1 to N
        for (i = i2 + 1; i <= n; i = i + 1) {
            tmp = a[(i - 1) + (i1 - 1) * lda];
            a[(i - 1) + (i1 - 1) * lda] = a[(i - 1) + (i2 - 1) * lda];
            a[(i - 1) + (i2 - 1) * lda] = tmp;
        }
        //
    }
    //
}
