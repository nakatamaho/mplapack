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

REAL RCsum1(INTEGER const n, COMPLEX *cx, INTEGER const incx) {
    REAL return_value = 0.0;
    REAL stemp = 0.0;
    INTEGER nincx = 0;
    INTEGER i = 0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    return_value = 0.0;
    stemp = 0.0;
    if (n <= 0) {
        return return_value;
    }
    if (incx == 1) {
        goto statement_20;
    }
    //
    //     CODE FOR INCREMENT NOT EQUAL TO 1
    //
    nincx = n * incx;
    for (i = 1; i <= nincx; i = i + incx) {
        //
        //        NEXT LINE MODIFIED.
        //
        stemp += abs(cx[i - 1]);
    }
    return_value = stemp;
    return return_value;
//
//     CODE FOR INCREMENT EQUAL TO 1
//
statement_20:
    for (i = 1; i <= n; i = i + 1) {
        //
        //        NEXT LINE MODIFIED.
        //
        stemp += abs(cx[i - 1]);
    }
    return_value = stemp;
    return return_value;
    //
    //     End of RCsum1
    //
}
