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

void Rlaord(const char *job, INTEGER const n, REAL *x, INTEGER const incx) {
    INTEGER inc = 0;
    INTEGER i = 0;
    INTEGER ix = 0;
    INTEGER ixnext = 0;
    REAL temp = 0.0;
    //
    //  -- LAPACK test routine --
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    inc = abs(incx);
    if (Mlsame(job, "I")) {
        //
        //        Sort in increasing order
        //
        for (i = 2; i <= n; i = i + 1) {
            ix = 1 + (i - 1) * inc;
        statement_10:
            if (ix == 1) {
                goto statement_20;
            }
            ixnext = ix - inc;
            if (x[ix - 1] > x[ixnext - 1]) {
                goto statement_20;
            } else {
                temp = x[ix - 1];
                x[ix - 1] = x[ixnext - 1];
                x[ixnext - 1] = temp;
            }
            ix = ixnext;
            goto statement_10;
        statement_20:;
        }
        //
    } else if (Mlsame(job, "D")) {
        //
        //        Sort in decreasing order
        //
        for (i = 2; i <= n; i = i + 1) {
            ix = 1 + (i - 1) * inc;
        statement_30:
            if (ix == 1) {
                goto statement_40;
            }
            ixnext = ix - inc;
            if (x[ix - 1] < x[ixnext - 1]) {
                goto statement_40;
            } else {
                temp = x[ix - 1];
                x[ix - 1] = x[ixnext - 1];
                x[ixnext - 1] = temp;
            }
            ix = ixnext;
            goto statement_30;
        statement_40:;
        }
    }
    //
    //     End of Rlaord
    //
}
