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

inline REAL cabs1(COMPLEX zdum) { return abs(zdum.real()) + abs(zdum.imag()); }

void Cpocon(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL const anorm, REAL &rcond, COMPLEX *work, REAL *rwork, INTEGER &info) {
    COMPLEX zdum = 0.0;
    bool upper = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL smlnum = 0.0;
    INTEGER kase = 0;
    char normin;
    REAL ainvnm = 0.0;
    INTEGER isave[3];
    REAL scalel = 0.0;
    REAL scaleu = 0.0;
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (anorm < zero) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Cpocon", -info);
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
    smlnum = Rlamch("Safe minimum");
    //
    //     Estimate the 1-norm of inv(A).
    //
    kase = 0;
    normin = 'N';
statement_10:
    Clacn2(n, &work[(n + 1) - 1], work, ainvnm, kase, isave);
    if (kase != 0) {
        if (upper) {
            //
            //           Multiply by inv(U**H).
            //
            Clatrs("Upper", "Conjugate transpose", "Non-unit", &normin, n, a, lda, work, scalel, rwork, info);
            normin = 'Y';
            //
            //           Multiply by inv(U).
            //
            Clatrs("Upper", "No transpose", "Non-unit", &normin, n, a, lda, work, scaleu, rwork, info);
        } else {
            //
            //           Multiply by inv(L).
            //
            Clatrs("Lower", "No transpose", "Non-unit", &normin, n, a, lda, work, scalel, rwork, info);
            normin = 'Y';
            //
            //           Multiply by inv(L**H).
            //
            Clatrs("Lower", "Conjugate transpose", "Non-unit", &normin, n, a, lda, work, scaleu, rwork, info);
        }
        //
        //        Multiply by 1/SCALE if doing so will not cause overflow.
        //
        scale = scalel * scaleu;
        if (scale != one) {
            ix = iCamax(n, work, 1);
            if (scale < cabs1(work[ix - 1]) * smlnum || scale == zero) {
                goto statement_20;
            }
            CRrscl(n, scale, work, 1);
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
    //     End of Cpocon
    //
}
