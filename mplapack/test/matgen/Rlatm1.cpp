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

void Rlatm1(INTEGER const mode, REAL const cond, INTEGER const irsign, INTEGER const idist, INTEGER *iseed, REAL *d, INTEGER const n, INTEGER &info) {
    const REAL one = 1.0;
    INTEGER i = 0;
    REAL alpha = 0.0;
    REAL temp = 0.0;
    const REAL half = 0.5e0;
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
    //     Decode and Test the input parameters. Initialize flags & seed.
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Set INFO if an error
    //
    if (mode < -6 || mode > 6) {
        info = -1;
    } else if ((mode != -6 && mode != 0 && mode != 6) && (irsign != 0 && irsign != 1)) {
        info = -2;
    } else if ((mode != -6 && mode != 0 && mode != 6) && cond < one) {
        info = -3;
    } else if ((mode == 6 || mode == -6) && (idist < 1 || idist > 3)) {
        info = -4;
    } else if (n < 0) {
        info = -7;
    }
    //
    if (info != 0) {
        Mxerbla("Rlatm1", -info);
        return;
    }
    //
    //     Compute D according to COND and MODE
    //
    if (mode != 0) {
        switch (abs(mode)) {
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
        default:
            break;
        }
    //
    //        One large D value:
    //
    statement_10:
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = one / cond;
        }
        d[1 - 1] = one;
        goto statement_120;
    //
    //        One small D value:
    //
    statement_30:
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = one;
        }
        d[n - 1] = one / cond;
        goto statement_120;
    //
    //        Exponentially distributed D values:
    //
    statement_50:
        d[1 - 1] = one;
        if (n > 1) {
	    alpha = pow(cond, (-one / castREAL(n - 1)));
            for (i = 2; i <= n; i = i + 1) {
	      d[i - 1] = pow(alpha, (i - 1));
            }
        }
        goto statement_120;
    //
    //        Arithmetically distributed D values:
    //
    statement_70:
        d[1 - 1] = one;
        if (n > 1) {
            temp = one / cond;
            alpha = (one - temp) / castREAL(n - 1);
            for (i = 2; i <= n; i = i + 1) {
                d[i - 1] = castREAL(n - i) * alpha + temp;
            }
        }
        goto statement_120;
    //
    //        Randomly distributed D values on ( 1/COND , 1):
    //
    statement_90:
        alpha = log(one / cond);
        for (i = 1; i <= n; i = i + 1) {
	    d[i - 1] = exp(alpha * Rlaran(iseed);
        }
        goto statement_120;
    //
    //        Randomly distributed D values from IDIST
    //
    statement_110:
        Rlarnv(idist, iseed, n, d);
    //
    statement_120:
        //
        //        If MODE neither -6 nor 0 nor 6, and IRSIGN = 1, assign
        //        random signs to D
        //
        if ((mode != -6 && mode != 0 && mode != 6) && irsign == 1) {
            for (i = 1; i <= n; i = i + 1) {
                temp = Rlaran[iseed - 1];
                if (temp > half) {
                    d[i - 1] = -d[i - 1];
                }
            }
        }
        //
        //        Reverse if MODE < 0
        //
        if (mode < 0) {
            for (i = 1; i <= n / 2; i = i + 1) {
                temp = d[i - 1];
                d[i - 1] = d[(n + 1 - i) - 1];
                d[(n + 1 - i) - 1] = temp;
            }
        }
        //
    }
    //
    //     End of Rlatm1
    //
}
