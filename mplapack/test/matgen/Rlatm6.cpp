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

#include <mplapack_matgen.h>

void Rlatm6(INTEGER const type, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, REAL *x, INTEGER const ldx, REAL *y, INTEGER const ldy, REAL const alpha, REAL const beta, REAL const wx, REAL const wy, REAL *s, REAL *dif) {
    //
    //  -- LAPACK computational routine --
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
    //     .. Local Arrays ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Generate test problem ...
    //     (Da, Db) ...
    //
    INTEGER i = 0;
    INTEGER j = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER ldb = lda;
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            //
            if (i == j) {
                a[(i - 1) + (i - 1) * lda] = castREAL(i) + alpha;
                b[(i - 1) + (i - 1) * ldb] = one;
            } else {
                a[(i - 1) + (j - 1) * lda] = zero;
                b[(i - 1) + (j - 1) * ldb] = zero;
            }
            //
        }
    }
    //
    //     Form X and Y
    //
    Rlacpy("F", n, n, b, lda, y, ldy);
    y[(3 - 1)] = -wy;
    y[(4 - 1)] = wy;
    y[(5 - 1)] = -wy;
    y[(3 - 1) + (2 - 1) * ldy] = -wy;
    y[(4 - 1) + (2 - 1) * ldy] = wy;
    y[(5 - 1) + (2 - 1) * ldy] = -wy;
    //
    Rlacpy("F", n, n, b, lda, x, ldx);
    x[(3 - 1) * ldx] = -wx;
    x[(4 - 1) * ldx] = -wx;
    x[(5 - 1) * ldx] = wx;
    x[(2 - 1) + (3 - 1) * ldx] = wx;
    x[(2 - 1) + (4 - 1) * ldx] = -wx;
    x[(2 - 1) + (5 - 1) * ldx] = -wx;
    //
    //     Form (A, B)
    //
    b[(3 - 1) * ldb] = wx + wy;
    b[(2 - 1) + (3 - 1) * ldb] = -wx + wy;
    b[(4 - 1) * ldb] = wx - wy;
    b[(2 - 1) + (4 - 1) * ldb] = wx - wy;
    b[(5 - 1) * ldb] = -wx + wy;
    b[(2 - 1) + (5 - 1) * ldb] = wx + wy;
    const REAL two = 2.0e+0;
    if (type == 1) {
        a[(3 - 1) * lda] = wx * a[(1 - 1)] + wy * a[(3 - 1) + (3 - 1) * lda];
        a[(2 - 1) + (3 - 1) * lda] = -wx * a[(2 - 1) + (2 - 1) * lda] + wy * a[(3 - 1) + (3 - 1) * lda];
        a[(4 - 1) * lda] = wx * a[(1 - 1)] - wy * a[(4 - 1) + (4 - 1) * lda];
        a[(2 - 1) + (4 - 1) * lda] = wx * a[(2 - 1) + (2 - 1) * lda] - wy * a[(4 - 1) + (4 - 1) * lda];
        a[(5 - 1) * lda] = -wx * a[(1 - 1)] + wy * a[(5 - 1) + (5 - 1) * lda];
        a[(2 - 1) + (5 - 1) * lda] = wx * a[(2 - 1) + (2 - 1) * lda] + wy * a[(5 - 1) + (5 - 1) * lda];
    } else if (type == 2) {
        a[(3 - 1) * lda] = two * wx + wy;
        a[(2 - 1) + (3 - 1) * lda] = wy;
        a[(4 - 1) * lda] = -wy * (two + alpha + beta);
        a[(2 - 1) + (4 - 1) * lda] = two * wx - wy * (two + alpha + beta);
        a[(5 - 1) * lda] = -two * wx + wy * (alpha - beta);
        a[(2 - 1) + (5 - 1) * lda] = wy * (alpha - beta);
        a[(1 - 1)] = one;
        a[(2 - 1) * lda] = -one;
        a[(2 - 1)] = one;
        a[(2 - 1) + (2 - 1) * lda] = a[(1 - 1)];
        a[(3 - 1) + (3 - 1) * lda] = one;
        a[(4 - 1) + (4 - 1) * lda] = one + alpha;
        a[(4 - 1) + (5 - 1) * lda] = one + beta;
        a[(5 - 1) + (4 - 1) * lda] = -a[(4 - 1) + (5 - 1) * lda];
        a[(5 - 1) + (5 - 1) * lda] = a[(4 - 1) + (4 - 1) * lda];
    }
    //
    //     Compute condition numbers
    //
    const REAL three = 3.0e+0;
    REAL z[12 * 12];
    REAL work[100];
    INTEGER info = 0;
    if (type == 1) {
        //
        s[1 - 1] = one / sqrt((one + three * wy * wy) / (one + a[(1 - 1)] * a[(1 - 1)]));
        s[2 - 1] = one / sqrt((one + three * wy * wy) / (one + a[(2 - 1) + (2 - 1) * lda] * a[(2 - 1) + (2 - 1) * lda]));
        s[3 - 1] = one / sqrt((one + two * wx * wx) / (one + a[(3 - 1) + (3 - 1) * lda] * a[(3 - 1) + (3 - 1) * lda]));
        s[4 - 1] = one / sqrt((one + two * wx * wx) / (one + a[(4 - 1) + (4 - 1) * lda] * a[(4 - 1) + (4 - 1) * lda]));
        s[5 - 1] = one / sqrt((one + two * wx * wx) / (one + a[(5 - 1) + (5 - 1) * lda] * a[(5 - 1) + (5 - 1) * lda]));
        //
        Rlakf2(1, 4, a, lda, &a[(2 - 1) + (2 - 1) * lda], b, &b[(2 - 1) + (2 - 1) * ldb], z, 12);
        Rgesvd("N", "N", 8, 8, z, 12, work, &work[9 - 1], 1, &work[10 - 1], 1, &work[11 - 1], 40, info);
        dif[1 - 1] = work[8 - 1];
        //
        Rlakf2(4, 1, a, lda, &a[(5 - 1) + (5 - 1) * lda], b, &b[(5 - 1) + (5 - 1) * ldb], z, 12);
        Rgesvd("N", "N", 8, 8, z, 12, work, &work[9 - 1], 1, &work[10 - 1], 1, &work[11 - 1], 40, info);
        dif[5 - 1] = work[8 - 1];
        //
    } else if (type == 2) {
        //
        s[1 - 1] = one / sqrt(one / three + wy * wy);
        s[2 - 1] = s[1 - 1];
        s[3 - 1] = one / sqrt(one / two + wx * wx);
        s[4 - 1] = one / sqrt((one + two * wx * wx) / (one + (one + alpha) * (one + alpha) + (one + beta) * (one + beta)));
        s[5 - 1] = s[4 - 1];
        //
        Rlakf2(2, 3, a, lda, &a[(3 - 1) + (3 - 1) * lda], b, &b[(3 - 1) + (3 - 1) * ldb], z, 12);
        Rgesvd("N", "N", 12, 12, z, 12, work, &work[13 - 1], 1, &work[14 - 1], 1, &work[15 - 1], 60, info);
        dif[1 - 1] = work[12 - 1];
        //
        Rlakf2(3, 2, a, lda, &a[(4 - 1) + (4 - 1) * lda], b, &b[(4 - 1) + (4 - 1) * ldb], z, 12);
        Rgesvd("N", "N", 12, 12, z, 12, work, &work[13 - 1], 1, &work[14 - 1], 1, &work[15 - 1], 60, info);
        dif[5 - 1] = work[12 - 1];
        //
    }
    //
    //     End of Rlatm6
    //
}
