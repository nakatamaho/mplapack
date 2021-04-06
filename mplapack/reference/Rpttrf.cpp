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

void Rpttrf(INTEGER const n, REAL *d, REAL *e, INTEGER &info) {
    INTEGER i4 = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL ei = 0.0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    if (n < 0) {
        info = -1;
        Mxerbla("Rpttrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Compute the L*D*L**T (or U**T*D*U) factorization of A.
    //
    i4 = mod(n - 1, 4);
    for (i = 1; i <= i4; i = i + 1) {
        if (d[i - 1] <= zero) {
            info = i;
            goto statement_30;
        }
        ei = e[i - 1];
        e[i - 1] = ei / d[i - 1];
        d[(i + 1) - 1] = d[(i + 1) - 1] - e[i - 1] * ei;
    }
    //
    for (i = i4 + 1; i <= n - 4; i = i + 4) {
        //
        //        Drop out of the loop if d(i) <= 0: the matrix is not positive
        //        definite.
        //
        if (d[i - 1] <= zero) {
            info = i;
            goto statement_30;
        }
        //
        //        Solve for e(i) and d(i+1).
        //
        ei = e[i - 1];
        e[i - 1] = ei / d[i - 1];
        d[(i + 1) - 1] = d[(i + 1) - 1] - e[i - 1] * ei;
        //
        if (d[(i + 1) - 1] <= zero) {
            info = i + 1;
            goto statement_30;
        }
        //
        //        Solve for e(i+1) and d(i+2).
        //
        ei = e[(i + 1) - 1];
        e[(i + 1) - 1] = ei / d[(i + 1) - 1];
        d[(i + 2) - 1] = d[(i + 2) - 1] - e[(i + 1) - 1] * ei;
        //
        if (d[(i + 2) - 1] <= zero) {
            info = i + 2;
            goto statement_30;
        }
        //
        //        Solve for e(i+2) and d(i+3).
        //
        ei = e[(i + 2) - 1];
        e[(i + 2) - 1] = ei / d[(i + 2) - 1];
        d[(i + 3) - 1] = d[(i + 3) - 1] - e[(i + 2) - 1] * ei;
        //
        if (d[(i + 3) - 1] <= zero) {
            info = i + 3;
            goto statement_30;
        }
        //
        //        Solve for e(i+3) and d(i+4).
        //
        ei = e[(i + 3) - 1];
        e[(i + 3) - 1] = ei / d[(i + 3) - 1];
        d[(i + 4) - 1] = d[(i + 4) - 1] - e[(i + 3) - 1] * ei;
    }
    //
    //     Check d(n) for positive definiteness.
    //
    if (d[n - 1] <= zero) {
        info = n;
    }
//
statement_30:;
    //
    //     End of Rpttrf
    //
}
