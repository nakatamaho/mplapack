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

void Rrscl(INTEGER const &n, REAL const &sa, REAL *sx, INTEGER const &incx) {
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL cden = 0.0;
    REAL cnum = 0.0;
    REAL cden1 = 0.0;
    REAL cnum1 = 0.0;
    const REAL zero = 0.0;
    REAL mul = 0.0;
    bool done = false;
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
    // =====================================================================
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
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     Get machine parameters
    //
    smlnum = dlamch("S");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Initialize the denominator to SA and the numerator to 1.
    //
    cden = sa;
    cnum = one;
//
statement_10:
    cden1 = cden * smlnum;
    cnum1 = cnum / bignum;
    if (abs(cden1) > abs(cnum) && cnum != zero) {
        //
        //        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
        //
        mul = smlnum;
        done = false;
        cden = cden1;
    } else if (abs(cnum1) > abs(cden)) {
        //
        //        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
        //
        mul = bignum;
        done = false;
        cnum = cnum1;
    } else {
        //
        //        Multiply X by CNUM / CDEN and return.
        //
        mul = cnum / cden;
        done = true;
    }
    //
    //     Scale the vector X by MUL
    //
    Rscal(n, mul, sx, incx);
    //
    if (!done) {
        goto statement_10;
    }
    //
    //     End of Rrscl
    //
}
