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

void Rtpt06(REAL const rcond, REAL const rcondc, const char *uplo, const char *diag, INTEGER const n, REAL *ap, REAL *work, REAL &rat) {
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL eps = Rlamch("Epsilon");
    REAL rmax = max(rcond, rcondc);
    REAL rmin = min(rcond, rcondc);
    //
    //     Do the easy cases first.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    REAL anorm = 0.0;
    if (rmin < zero) {
        //
        //        Invalid value for RCOND or RCONDC, return 1/EPS.
        //
        rat = one / eps;
        //
    } else if (rmin > zero) {
        //
        //        Both estimates are positive, return RMAX/RMIN - 1.
        //
        rat = rmax / rmin - one;
        //
    } else if (rmax == zero) {
        //
        //        Both estimates zero.
        //
        rat = zero;
        //
    } else {
        //
        //        One estimate is zero, the other is non-zero.  If the matrix is
        //        ill-conditioned, return the nonzero estimate multiplied by
        //        1/EPS; if the matrix is badly scaled, return the nonzero
        //        estimate multiplied by BIGNUM/TMAX, where TMAX is the maximum
        //        element in absolute value in A.
        //
        smlnum = Rlamch("Safe minimum");
        bignum = one / smlnum;
        Rlabad(smlnum, bignum);
        anorm = Rlantp("M", uplo, diag, n, ap, work);
        //
        rat = rmax * (min(bignum / max(one, anorm), one / eps));
    }
    //
    //     End of Rtpt06
    //
}
