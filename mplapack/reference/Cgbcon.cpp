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

inline REAL abs1(COMPLEX zdum) { return (abs(zdum.real()) + abs(zdum.imag())); }

void Cgbcon(const char *norm, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX *ab, INTEGER const ldab, INTEGER *ipiv, REAL const anorm, REAL &rcond, COMPLEX *work, REAL *rwork, INTEGER &info) {
    COMPLEX zdum = 0.0;
    bool onenrm = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL smlnum = 0.0;
    REAL ainvnm = 0.0;
    char normin ;
    INTEGER kase1 = 0;
    INTEGER kd = 0;
    bool lnoti = false;
    INTEGER kase = 0;
    INTEGER isave[3];
    INTEGER j = 0;
    INTEGER lm = 0;
    INTEGER jp = 0;
    COMPLEX t = 0.0;
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
    onenrm = Mlsame(norm , "1") || Mlsame(norm, "O");
    if (!onenrm && !Mlsame(norm, "I")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0) {
        info = -3;
    } else if (ku < 0) {
        info = -4;
    } else if (ldab < 2 * kl + ku + 1) {
        info = -6;
    } else if (anorm < zero) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Cgbcon", -info);
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
    //     Estimate the norm of inv(A).
    //
    ainvnm = zero;
    normin = 'N';
    if (onenrm) {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kd = kl + ku + 1;
    lnoti = kl > 0;
    kase = 0;
statement_10:
    Clacn2(n, &work[(n + 1) - 1], work, ainvnm, kase, isave);
    if (kase != 0) {
        if (kase == kase1) {
            //
            //           Multiply by inv(L).
            //
            if (lnoti) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    lm = min(kl, n - j);
                    jp = ipiv[j - 1];
                    t = work[jp - 1];
                    if (jp != j) {
                        work[jp - 1] = work[j - 1];
                        work[j - 1] = t;
                    }
                    Caxpy(lm, -t, &ab[((kd + 1) - 1) + (j - 1) * ldab], 1, &work[(j + 1) - 1], 1);
                }
            }
            //
            //           Multiply by inv(U).
            //
            Clatbs("Upper", "No transpose", "Non-unit", &normin, n, kl + ku, ab, ldab, work, scale, rwork, info);
        } else {
            //
            //           Multiply by inv(U**H).
            //
            Clatbs("Upper", "Conjugate transpose", "Non-unit", &normin, n, kl + ku, ab, ldab, work, scale, rwork, info);
            //
            //           Multiply by inv(L**H).
            //
            if (lnoti) {
                for (j = n - 1; j >= 1; j = j - 1) {
                    lm = min(kl, n - j);
                    work[j - 1] = work[j - 1] - Cdotc(lm, &ab[((kd + 1) - 1) + (j - 1) * ldab], 1, &work[(j + 1) - 1], 1);
                    jp = ipiv[j - 1];
                    if (jp != j) {
                        t = work[jp - 1];
                        work[jp - 1] = work[j - 1];
                        work[j - 1] = t;
                    }
                }
            }
        }
        //
        //        Divide X by 1/SCALE if doing so will not cause overflow.
        //
        normin = 'Y';
        if (scale != one) {
            ix = iCamax(n, work, 1);
            if (scale < abs1(work[ix - 1]) * smlnum || scale == zero) {
                goto statement_40;
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
statement_40:;
    //
    //     End of Cgbcon
    //
}
