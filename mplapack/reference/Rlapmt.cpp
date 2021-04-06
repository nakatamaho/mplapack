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

void Rlapmt(bool const forwrd, INTEGER const m, INTEGER const n, REAL *x, INTEGER const ldx, INTEGER *k) {
    INTEGER i = 0;
    INTEGER j = 0;
    INTEGER in = 0;
    INTEGER ii = 0;
    REAL temp = 0.0;
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
    //     .. Local Scalars ..
    //     ..
    //     .. Executable Statements ..
    //
    if (n <= 1) {
        return;
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        k[i - 1] = -k[i - 1];
    }
    //
    if (forwrd) {
        //
        //        Forward permutation
        //
        for (i = 1; i <= n; i = i + 1) {
            //
            if (k[i - 1] > 0) {
                goto statement_40;
            }
            //
            j = i;
            k[j - 1] = -k[j - 1];
            in = k[j - 1];
        //
        statement_20:
            if (k[in - 1] > 0) {
                goto statement_40;
            }
            //
            for (ii = 1; ii <= m; ii = ii + 1) {
                temp = x[(ii - 1) + (j - 1) * ldx];
                x[(ii - 1) + (j - 1) * ldx] = x[(ii - 1) + (in - 1) * ldx];
                x[(ii - 1) + (in - 1) * ldx] = temp;
            }
            //
            k[in - 1] = -k[in - 1];
            j = in;
            in = k[in - 1];
            goto statement_20;
        //
        statement_40:;
            //
        }
        //
    } else {
        //
        //        Backward permutation
        //
        for (i = 1; i <= n; i = i + 1) {
            //
            if (k[i - 1] > 0) {
                goto statement_80;
            }
            //
            k[i - 1] = -k[i - 1];
            j = k[i - 1];
        statement_60:
            if (j == i) {
                goto statement_80;
            }
            //
            for (ii = 1; ii <= m; ii = ii + 1) {
                temp = x[(ii - 1) + (i - 1) * ldx];
                x[(ii - 1) + (i - 1) * ldx] = x[(ii - 1) + (j - 1) * ldx];
                x[(ii - 1) + (j - 1) * ldx] = temp;
            }
            //
            k[j - 1] = -k[j - 1];
            j = k[j - 1];
            goto statement_60;
        //
        statement_80:;
            //
        }
        //
    }
    //
    //     End of Rlapmt
    //
}
