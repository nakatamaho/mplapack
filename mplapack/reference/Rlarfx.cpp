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

void Rlarfx(const char *side, INTEGER const m, INTEGER const n, REAL *v, REAL const tau, REAL *c, INTEGER const ldc, REAL *work) {
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL t1 = 0.0;
    INTEGER j = 0;
    REAL v1 = 0.0;
    REAL v2 = 0.0;
    REAL t2 = 0.0;
    REAL sum = 0.0;
    REAL v3 = 0.0;
    REAL t3 = 0.0;
    REAL v4 = 0.0;
    REAL t4 = 0.0;
    REAL v5 = 0.0;
    REAL t5 = 0.0;
    REAL v6 = 0.0;
    REAL t6 = 0.0;
    REAL v7 = 0.0;
    REAL t7 = 0.0;
    REAL v8 = 0.0;
    REAL t8 = 0.0;
    REAL v9 = 0.0;
    REAL t9 = 0.0;
    REAL v10 = 0.0;
    REAL t10 = 0.0;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    if (tau == zero) {
        return;
    }
    if (Mlsame(side, "L")) {
        //
        //        Form  H * C, where H has order m.
        //
        switch (m) {
        case 1:
            goto statement_10;
        case 2:
            goto statement_30;
        case 3:
            goto statement_50;
        case 4:
            goto statement_70;
        case 5:
            goto statement_90;
        case 6:
            goto statement_110;
        case 7:
            goto statement_130;
        case 8:
            goto statement_150;
        case 9:
            goto statement_170;
        case 10:
            goto statement_190;
        default:
            break;
        }
        //
        //        Code for general M
        //
        Rlarf(side, m, n, v, 1, tau, c, ldc, work);
        goto statement_410;
    statement_10:
        //
        //        Special code for 1 x 1 Householder
        //
        t1 = one - tau * v[1 - 1] * v[1 - 1];
        for (j = 1; j <= n; j = j + 1) {
            c[(j - 1) * ldc] = t1 * c[(j - 1) * ldc];
        }
        goto statement_410;
    statement_30:
        //
        //        Special code for 2 x 2 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
        }
        goto statement_410;
    statement_50:
        //
        //        Special code for 3 x 3 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
        }
        goto statement_410;
    statement_70:
        //
        //        Special code for 4 x 4 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
        }
        goto statement_410;
    statement_90:
        //
        //        Special code for 5 x 5 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc] + v5 * c[(5 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
            c[(5 - 1) + (j - 1) * ldc] = c[(5 - 1) + (j - 1) * ldc] - sum * t5;
        }
        goto statement_410;
    statement_110:
        //
        //        Special code for 6 x 6 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc] + v5 * c[(5 - 1) + (j - 1) * ldc] + v6 * c[(6 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
            c[(5 - 1) + (j - 1) * ldc] = c[(5 - 1) + (j - 1) * ldc] - sum * t5;
            c[(6 - 1) + (j - 1) * ldc] = c[(6 - 1) + (j - 1) * ldc] - sum * t6;
        }
        goto statement_410;
    statement_130:
        //
        //        Special code for 7 x 7 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc] + v5 * c[(5 - 1) + (j - 1) * ldc] + v6 * c[(6 - 1) + (j - 1) * ldc] + v7 * c[(7 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
            c[(5 - 1) + (j - 1) * ldc] = c[(5 - 1) + (j - 1) * ldc] - sum * t5;
            c[(6 - 1) + (j - 1) * ldc] = c[(6 - 1) + (j - 1) * ldc] - sum * t6;
            c[(7 - 1) + (j - 1) * ldc] = c[(7 - 1) + (j - 1) * ldc] - sum * t7;
        }
        goto statement_410;
    statement_150:
        //
        //        Special code for 8 x 8 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        v8 = v[8 - 1];
        t8 = tau * v8;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc] + v5 * c[(5 - 1) + (j - 1) * ldc] + v6 * c[(6 - 1) + (j - 1) * ldc] + v7 * c[(7 - 1) + (j - 1) * ldc] + v8 * c[(8 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
            c[(5 - 1) + (j - 1) * ldc] = c[(5 - 1) + (j - 1) * ldc] - sum * t5;
            c[(6 - 1) + (j - 1) * ldc] = c[(6 - 1) + (j - 1) * ldc] - sum * t6;
            c[(7 - 1) + (j - 1) * ldc] = c[(7 - 1) + (j - 1) * ldc] - sum * t7;
            c[(8 - 1) + (j - 1) * ldc] = c[(8 - 1) + (j - 1) * ldc] - sum * t8;
        }
        goto statement_410;
    statement_170:
        //
        //        Special code for 9 x 9 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        v8 = v[8 - 1];
        t8 = tau * v8;
        v9 = v[9 - 1];
        t9 = tau * v9;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc] + v5 * c[(5 - 1) + (j - 1) * ldc] + v6 * c[(6 - 1) + (j - 1) * ldc] + v7 * c[(7 - 1) + (j - 1) * ldc] + v8 * c[(8 - 1) + (j - 1) * ldc] + v9 * c[(9 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
            c[(5 - 1) + (j - 1) * ldc] = c[(5 - 1) + (j - 1) * ldc] - sum * t5;
            c[(6 - 1) + (j - 1) * ldc] = c[(6 - 1) + (j - 1) * ldc] - sum * t6;
            c[(7 - 1) + (j - 1) * ldc] = c[(7 - 1) + (j - 1) * ldc] - sum * t7;
            c[(8 - 1) + (j - 1) * ldc] = c[(8 - 1) + (j - 1) * ldc] - sum * t8;
            c[(9 - 1) + (j - 1) * ldc] = c[(9 - 1) + (j - 1) * ldc] - sum * t9;
        }
        goto statement_410;
    statement_190:
        //
        //        Special code for 10 x 10 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        v8 = v[8 - 1];
        t8 = tau * v8;
        v9 = v[9 - 1];
        t9 = tau * v9;
        v10 = v[10 - 1];
        t10 = tau * v10;
        for (j = 1; j <= n; j = j + 1) {
            sum = v1 * c[(j - 1) * ldc] + v2 * c[(2 - 1) + (j - 1) * ldc] + v3 * c[(3 - 1) + (j - 1) * ldc] + v4 * c[(4 - 1) + (j - 1) * ldc] + v5 * c[(5 - 1) + (j - 1) * ldc] + v6 * c[(6 - 1) + (j - 1) * ldc] + v7 * c[(7 - 1) + (j - 1) * ldc] + v8 * c[(8 - 1) + (j - 1) * ldc] + v9 * c[(9 - 1) + (j - 1) * ldc] + v10 * c[(10 - 1) + (j - 1) * ldc];
            c[(j - 1) * ldc] = c[(j - 1) * ldc] - sum * t1;
            c[(2 - 1) + (j - 1) * ldc] = c[(2 - 1) + (j - 1) * ldc] - sum * t2;
            c[(3 - 1) + (j - 1) * ldc] = c[(3 - 1) + (j - 1) * ldc] - sum * t3;
            c[(4 - 1) + (j - 1) * ldc] = c[(4 - 1) + (j - 1) * ldc] - sum * t4;
            c[(5 - 1) + (j - 1) * ldc] = c[(5 - 1) + (j - 1) * ldc] - sum * t5;
            c[(6 - 1) + (j - 1) * ldc] = c[(6 - 1) + (j - 1) * ldc] - sum * t6;
            c[(7 - 1) + (j - 1) * ldc] = c[(7 - 1) + (j - 1) * ldc] - sum * t7;
            c[(8 - 1) + (j - 1) * ldc] = c[(8 - 1) + (j - 1) * ldc] - sum * t8;
            c[(9 - 1) + (j - 1) * ldc] = c[(9 - 1) + (j - 1) * ldc] - sum * t9;
            c[(10 - 1) + (j - 1) * ldc] = c[(10 - 1) + (j - 1) * ldc] - sum * t10;
        }
        goto statement_410;
    } else {
        //
        //        Form  C * H, where H has order n.
        //
        switch (n) {
        case 1:
            goto statement_210;
        case 2:
            goto statement_230;
        case 3:
            goto statement_250;
        case 4:
            goto statement_270;
        case 5:
            goto statement_290;
        case 6:
            goto statement_310;
        case 7:
            goto statement_330;
        case 8:
            goto statement_350;
        case 9:
            goto statement_370;
        case 10:
            goto statement_390;
        default:
            break;
        }
        //
        //        Code for general N
        //
        Rlarf(side, m, n, v, 1, tau, c, ldc, work);
        goto statement_410;
    statement_210:
        //
        //        Special code for 1 x 1 Householder
        //
        t1 = one - tau * v[1 - 1] * v[1 - 1];
        for (j = 1; j <= m; j = j + 1) {
            c[(j - 1)] = t1 * c[(j - 1)];
        }
        goto statement_410;
    statement_230:
        //
        //        Special code for 2 x 2 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
        }
        goto statement_410;
    statement_250:
        //
        //        Special code for 3 x 3 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
        }
        goto statement_410;
    statement_270:
        //
        //        Special code for 4 x 4 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
        }
        goto statement_410;
    statement_290:
        //
        //        Special code for 5 x 5 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc] + v5 * c[(j - 1) + (5 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
            c[(j - 1) + (5 - 1) * ldc] = c[(j - 1) + (5 - 1) * ldc] - sum * t5;
        }
        goto statement_410;
    statement_310:
        //
        //        Special code for 6 x 6 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc] + v5 * c[(j - 1) + (5 - 1) * ldc] + v6 * c[(j - 1) + (6 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
            c[(j - 1) + (5 - 1) * ldc] = c[(j - 1) + (5 - 1) * ldc] - sum * t5;
            c[(j - 1) + (6 - 1) * ldc] = c[(j - 1) + (6 - 1) * ldc] - sum * t6;
        }
        goto statement_410;
    statement_330:
        //
        //        Special code for 7 x 7 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc] + v5 * c[(j - 1) + (5 - 1) * ldc] + v6 * c[(j - 1) + (6 - 1) * ldc] + v7 * c[(j - 1) + (7 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
            c[(j - 1) + (5 - 1) * ldc] = c[(j - 1) + (5 - 1) * ldc] - sum * t5;
            c[(j - 1) + (6 - 1) * ldc] = c[(j - 1) + (6 - 1) * ldc] - sum * t6;
            c[(j - 1) + (7 - 1) * ldc] = c[(j - 1) + (7 - 1) * ldc] - sum * t7;
        }
        goto statement_410;
    statement_350:
        //
        //        Special code for 8 x 8 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        v8 = v[8 - 1];
        t8 = tau * v8;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc] + v5 * c[(j - 1) + (5 - 1) * ldc] + v6 * c[(j - 1) + (6 - 1) * ldc] + v7 * c[(j - 1) + (7 - 1) * ldc] + v8 * c[(j - 1) + (8 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
            c[(j - 1) + (5 - 1) * ldc] = c[(j - 1) + (5 - 1) * ldc] - sum * t5;
            c[(j - 1) + (6 - 1) * ldc] = c[(j - 1) + (6 - 1) * ldc] - sum * t6;
            c[(j - 1) + (7 - 1) * ldc] = c[(j - 1) + (7 - 1) * ldc] - sum * t7;
            c[(j - 1) + (8 - 1) * ldc] = c[(j - 1) + (8 - 1) * ldc] - sum * t8;
        }
        goto statement_410;
    statement_370:
        //
        //        Special code for 9 x 9 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        v8 = v[8 - 1];
        t8 = tau * v8;
        v9 = v[9 - 1];
        t9 = tau * v9;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc] + v5 * c[(j - 1) + (5 - 1) * ldc] + v6 * c[(j - 1) + (6 - 1) * ldc] + v7 * c[(j - 1) + (7 - 1) * ldc] + v8 * c[(j - 1) + (8 - 1) * ldc] + v9 * c[(j - 1) + (9 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
            c[(j - 1) + (5 - 1) * ldc] = c[(j - 1) + (5 - 1) * ldc] - sum * t5;
            c[(j - 1) + (6 - 1) * ldc] = c[(j - 1) + (6 - 1) * ldc] - sum * t6;
            c[(j - 1) + (7 - 1) * ldc] = c[(j - 1) + (7 - 1) * ldc] - sum * t7;
            c[(j - 1) + (8 - 1) * ldc] = c[(j - 1) + (8 - 1) * ldc] - sum * t8;
            c[(j - 1) + (9 - 1) * ldc] = c[(j - 1) + (9 - 1) * ldc] - sum * t9;
        }
        goto statement_410;
    statement_390:
        //
        //        Special code for 10 x 10 Householder
        //
        v1 = v[1 - 1];
        t1 = tau * v1;
        v2 = v[2 - 1];
        t2 = tau * v2;
        v3 = v[3 - 1];
        t3 = tau * v3;
        v4 = v[4 - 1];
        t4 = tau * v4;
        v5 = v[5 - 1];
        t5 = tau * v5;
        v6 = v[6 - 1];
        t6 = tau * v6;
        v7 = v[7 - 1];
        t7 = tau * v7;
        v8 = v[8 - 1];
        t8 = tau * v8;
        v9 = v[9 - 1];
        t9 = tau * v9;
        v10 = v[10 - 1];
        t10 = tau * v10;
        for (j = 1; j <= m; j = j + 1) {
            sum = v1 * c[(j - 1)] + v2 * c[(j - 1) + (2 - 1) * ldc] + v3 * c[(j - 1) + (3 - 1) * ldc] + v4 * c[(j - 1) + (4 - 1) * ldc] + v5 * c[(j - 1) + (5 - 1) * ldc] + v6 * c[(j - 1) + (6 - 1) * ldc] + v7 * c[(j - 1) + (7 - 1) * ldc] + v8 * c[(j - 1) + (8 - 1) * ldc] + v9 * c[(j - 1) + (9 - 1) * ldc] + v10 * c[(j - 1) + (10 - 1) * ldc];
            c[(j - 1)] = c[(j - 1)] - sum * t1;
            c[(j - 1) + (2 - 1) * ldc] = c[(j - 1) + (2 - 1) * ldc] - sum * t2;
            c[(j - 1) + (3 - 1) * ldc] = c[(j - 1) + (3 - 1) * ldc] - sum * t3;
            c[(j - 1) + (4 - 1) * ldc] = c[(j - 1) + (4 - 1) * ldc] - sum * t4;
            c[(j - 1) + (5 - 1) * ldc] = c[(j - 1) + (5 - 1) * ldc] - sum * t5;
            c[(j - 1) + (6 - 1) * ldc] = c[(j - 1) + (6 - 1) * ldc] - sum * t6;
            c[(j - 1) + (7 - 1) * ldc] = c[(j - 1) + (7 - 1) * ldc] - sum * t7;
            c[(j - 1) + (8 - 1) * ldc] = c[(j - 1) + (8 - 1) * ldc] - sum * t8;
            c[(j - 1) + (9 - 1) * ldc] = c[(j - 1) + (9 - 1) * ldc] - sum * t9;
            c[(j - 1) + (10 - 1) * ldc] = c[(j - 1) + (10 - 1) * ldc] - sum * t10;
        }
        goto statement_410;
    }
statement_410:;
    //
    //     End of Rlarfx
    //
}
