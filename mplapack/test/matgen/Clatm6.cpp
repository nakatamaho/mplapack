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

void Clatm6(INTEGER const type, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, COMPLEX *x, INTEGER const ldx, COMPLEX *y, INTEGER const ldy, COMPLEX const alpha, COMPLEX const beta, COMPLEX const wx, COMPLEX const wy, REAL *s, REAL *dif) {
    //
    //     Generate test problem ...
    //     (Da, Db) ...
    //
    INTEGER i = 0;
    INTEGER j = 0;
    INTEGER ldb = lda;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            //
            if (i == j) {
                a[(i - 1) + (i - 1) * lda] = COMPLEX(i) + alpha;
                b[(i - 1) + (i - 1) * ldb] = one;
            } else {
                a[(i - 1) + (j - 1) * lda] = zero;
                b[(i - 1) + (j - 1) * ldb] = zero;
            }
            //
        }
    }
    const REAL rone = 1.0;
    if (type == 2) {
        a[(1 - 1)] = COMPLEX(rone, rone);
        a[(2 - 1) + (2 - 1) * lda] = conj(a[(1 - 1)]);
        a[(3 - 1) + (3 - 1) * lda] = one;
        a[(4 - 1) + (4 - 1) * lda] = COMPLEX((one + alpha).real(), (one + beta).real());
        a[(5 - 1) + (5 - 1) * lda] = conj(a[(4 - 1) + (4 - 1) * lda]);
    }
    //
    //     Form X and Y
    //
    Clacpy("F", n, n, b, lda, y, ldy);
    y[(3 - 1)] = -conj(wy);
    y[(4 - 1)] = conj(wy);
    y[(5 - 1)] = -conj(wy);
    y[(3 - 1) + (2 - 1) * ldy] = -conj(wy);
    y[(4 - 1) + (2 - 1) * ldy] = conj(wy);
    y[(5 - 1) + (2 - 1) * ldy] = -conj(wy);
    //
    Clacpy("F", n, n, b, lda, x, ldx);
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
    a[(3 - 1) * lda] = wx * a[(1 - 1)] + wy * a[(3 - 1) + (3 - 1) * lda];
    a[(2 - 1) + (3 - 1) * lda] = -wx * a[(2 - 1) + (2 - 1) * lda] + wy * a[(3 - 1) + (3 - 1) * lda];
    a[(4 - 1) * lda] = wx * a[(1 - 1)] - wy * a[(4 - 1) + (4 - 1) * lda];
    a[(2 - 1) + (4 - 1) * lda] = wx * a[(2 - 1) + (2 - 1) * lda] - wy * a[(4 - 1) + (4 - 1) * lda];
    a[(5 - 1) * lda] = -wx * a[(1 - 1)] + wy * a[(5 - 1) + (5 - 1) * lda];
    a[(2 - 1) + (5 - 1) * lda] = wx * a[(2 - 1) + (2 - 1) * lda] + wy * a[(5 - 1) + (5 - 1) * lda];
    //
    //     Compute condition numbers
    //
    const REAL three = 3.0e+0;
    s[1 - 1] = rone / sqrt((rone + three * abs(wy) * abs(wy)) / (rone + abs(a[(1 - 1)]) * abs(a[(1 - 1)])));
    s[2 - 1] = rone / sqrt((rone + three * abs(wy) * abs(wy)) / (rone + abs(a[(2 - 1) + (2 - 1) * lda]) * abs(a[(2 - 1) + (2 - 1) * lda])));
    const REAL two = 2.0e+0;
    s[3 - 1] = rone / sqrt((rone + two * abs(wx) * abs(wx)) / (rone + abs(a[(3 - 1) + (3 - 1) * lda]) * abs(a[(3 - 1) + (3 - 1) * lda])));
    s[4 - 1] = rone / sqrt((rone + two * abs(wx) * abs(wx)) / (rone + abs(a[(4 - 1) + (4 - 1) * lda]) * abs(a[(4 - 1) + (4 - 1) * lda])));
    s[5 - 1] = rone / sqrt((rone + two * abs(wx) * abs(wx)) / (rone + abs(a[(5 - 1) + (5 - 1) * lda]) * abs(a[(5 - 1) + (5 - 1) * lda])));
    //
    COMPLEX z[8 * 8];
    Clakf2(1, 4, a, lda, &a[(2 - 1) + (2 - 1) * lda], b, &b[(2 - 1) + (2 - 1) * ldb], z, 8);
    REAL rwork[50];
    COMPLEX work[26];
    INTEGER info = 0;
    Cgesvd("N", "N", 8, 8, z, 8, rwork, work, 1, &work[2 - 1], 1, &work[3 - 1], 24, &rwork[9 - 1], info);
    dif[1 - 1] = rwork[8 - 1];
    //
    Clakf2(4, 1, a, lda, &a[(5 - 1) + (5 - 1) * lda], b, &b[(5 - 1) + (5 - 1) * ldb], z, 8);
    Cgesvd("N", "N", 8, 8, z, 8, rwork, work, 1, &work[2 - 1], 1, &work[3 - 1], 24, &rwork[9 - 1], info);
    dif[5 - 1] = rwork[8 - 1];
    //
    //     End of Clatm6
    //
}
