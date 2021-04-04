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

void Rgecon(const char *norm, INTEGER const &n, REAL *a, INTEGER const &lda, REAL const &anorm, REAL &rcond, REAL *work, INTEGER *iwork, INTEGER &info) {
    bool onenrm = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL smlnum = 0.0;
    REAL ainvnm = 0.0;
    str<1> normin = char0;
    INTEGER kase1 = 0;
    INTEGER kase = 0;
    arr_1d<3, INTEGER> isave(fill0);
    REAL sl = 0.0;
    REAL su = 0.0;
    REAL scale = 0.0;
    INTEGER ix = 0;
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
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
    onenrm = norm == "1" || Mlsame(norm, "O");
    if (!onenrm && !Mlsame(norm, "I")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (anorm < zero) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Rgecon", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    rcond = zero;
    if (n == 0) {
        rcond = one;
        return;
    } else if (anorm == zero) {
        return;
    }
    //
    smlnum = dlamch("Safe minimum");
    //
    //     Estimate the norm of inv(A).
    //
    ainvnm = zero;
    normin = "N";
    if (onenrm) {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kase = 0;
statement_10:
    Rlacn2(n, work[(n + 1) - 1], work, iwork, ainvnm, kase, isave);
    if (kase != 0) {
        if (kase == kase1) {
            //
            //           Multiply by inv(L).
            //
            Rlatrs("Lower", "No transpose", "Unit", normin, n, a, lda, work, sl, work[(2 * n + 1) - 1], info);
            //
            //           Multiply by inv(U).
            //
            Rlatrs("Upper", "No transpose", "Non-unit", normin, n, a, lda, work, su, work[(3 * n + 1) - 1], info);
        } else {
            //
            //           Multiply by inv(U**T).
            //
            Rlatrs("Upper", "Transpose", "Non-unit", normin, n, a, lda, work, su, work[(3 * n + 1) - 1], info);
            //
            //           Multiply by inv(L**T).
            //
            Rlatrs("Lower", "Transpose", "Unit", normin, n, a, lda, work, sl, work[(2 * n + 1) - 1], info);
        }
        //
        //        Divide X by 1/(SL*SU) if doing so will not cause overflow.
        //
        scale = sl * su;
        normin = "Y";
        if (scale != one) {
            ix = iRamax[(n - 1) + (work - 1) * ldiRamax];
            if (scale < abs(work[ix - 1]) * smlnum || scale == zero) {
                goto statement_20;
            }
            Rrscl(n, scale, work, 1);
        }
        goto statement_10;
    }
    //
    //     Compute the estimate of the reciprocal condition number.
    //
    if (ainvnm != zero) {
        rcond = (one / ainvnm) / anorm;
    }
//
statement_20:;
    //
    //     End of Rgecon
    //
}
