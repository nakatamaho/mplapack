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

void Rlacn2(INTEGER const n, REAL *v, REAL *x, INTEGER *isgn, REAL &est, INTEGER &kase, INTEGER *isave) {
    INTEGER i = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    REAL estold = 0.0;
    REAL xs = 0.0;
    INTEGER jlast = 0;
    const INTEGER itmax = 5;
    REAL altsgn = 0.0;
    const REAL two = 2.0e+0;
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
    if (kase == 0) {
        for (i = 1; i <= n; i = i + 1) {
            x[i - 1] = one / n.real();
        }
        kase = 1;
        isave[1 - 1] = 1;
        return;
    }
    //
    switch (isave[1 - 1]) {
    case 1:
        goto statement_20;
    case 2:
        goto statement_40;
    case 3:
        goto statement_70;
    case 4:
        goto statement_110;
    case 5:
        goto statement_140;
    default:
        break;
    }
//
//     ................ ENTRY   (ISAVE( 1 ) = 1)
//     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
//
statement_20:
    if (n == 1) {
        v[1 - 1] = x[1 - 1];
        est = abs(v[1 - 1]);
        //        ... QUIT
        goto statement_150;
    }
    est = Rasum(n, x, 1);
    //
    for (i = 1; i <= n; i = i + 1) {
        if (x[i - 1] >= zero) {
            x[i - 1] = one;
        } else {
            x[i - 1] = -one;
        }
        isgn[i - 1] = nint[x[i - 1] - 1];
    }
    kase = 2;
    isave[1 - 1] = 2;
    return;
//
//     ................ ENTRY   (ISAVE( 1 ) = 2)
//     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
//
statement_40:
    isave[2 - 1] = iRamax(n, x, 1);
    isave[3 - 1] = 2;
//
//     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
//
statement_50:
    for (i = 1; i <= n; i = i + 1) {
        x[i - 1] = zero;
    }
    x[isave[2 - 1] - 1] = one;
    kase = 1;
    isave[1 - 1] = 3;
    return;
//
//     ................ ENTRY   (ISAVE( 1 ) = 3)
//     X HAS BEEN OVERWRITTEN BY A*X.
//
statement_70:
    Rcopy(n, x, 1, v, 1);
    estold = est;
    est = Rasum(n, v, 1);
    for (i = 1; i <= n; i = i + 1) {
        if (x[i - 1] >= zero) {
            xs = one;
        } else {
            xs = -one;
        }
        if (nint[xs - 1] != isgn[i - 1]) {
            goto statement_90;
        }
    }
    //     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
    goto statement_120;
//
statement_90:
    //     TEST FOR CYCLING.
    if (est <= estold) {
        goto statement_120;
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        if (x[i - 1] >= zero) {
            x[i - 1] = one;
        } else {
            x[i - 1] = -one;
        }
        isgn[i - 1] = nint[x[i - 1] - 1];
    }
    kase = 2;
    isave[1 - 1] = 4;
    return;
//
//     ................ ENTRY   (ISAVE( 1 ) = 4)
//     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
//
statement_110:
    jlast = isave[2 - 1];
    isave[2 - 1] = iRamax(n, x, 1);
    if ((x[jlast - 1] != abs(x[isave[2 - 1] - 1])) && (isave[3 - 1] < itmax)) {
        isave[3 - 1]++;
        goto statement_50;
    }
//
//     ITERATION COMPLETE.  FINAL STAGE.
//
statement_120:
    altsgn = one;
    for (i = 1; i <= n; i = i + 1) {
        x[i - 1] = altsgn * (one + i - 1.real() / n - 1.real());
        altsgn = -altsgn;
    }
    kase = 1;
    isave[1 - 1] = 5;
    return;
//
//     ................ ENTRY   (ISAVE( 1 ) = 5)
//     X HAS BEEN OVERWRITTEN BY A*X.
//
statement_140:
    temp = two * (Rasum(n, x, 1) / 3 * n.real());
    if (temp > est) {
        Rcopy(n, x, 1, v, 1);
        est = temp;
    }
//
statement_150:
    kase = 0;
    //
    //     End of Rlacn2
    //
}
