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

inline REAL abs1(COMPLEX zdum) { return abs(zdum.real()) + abs(zdum.imag()); }

REAL Cla_porcond_x(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *af, INTEGER const ldaf, COMPLEX *x, INTEGER &info, COMPLEX *work, REAL *rwork) {
    REAL return_value = 0.0;
    COMPLEX zdum = 0.0;
    bool upper = false;
    bool up = false;
    REAL anorm = 0.0;
    INTEGER i = 0;
    REAL tmp = 0.0;
    INTEGER j = 0;
    REAL ainvnm = 0.0;
    INTEGER kase = 0;
    INTEGER isave[3];
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
    //     .. Statement Function Definitions ..
    //     ..
    //     .. Executable Statements ..
    //
    return_value = 0.0;
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (ldaf < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Cla_porcond_x", -info);
        return return_value;
    }
    up = false;
    if (Mlsame(uplo, "U")) {
        up = true;
    }
    //
    //     Compute norm of op(A)*op2(C).
    //
    anorm = 0.0;
    if (up) {
        for (i = 1; i <= n; i = i + 1) {
            tmp = 0.0;
            for (j = 1; j <= i; j = j + 1) {
                tmp += abs1(a[(j - 1) + (i - 1) * lda] * x[j - 1]);
            }
            for (j = i + 1; j <= n; j = j + 1) {
                tmp += abs1(a[(i - 1) + (j - 1) * lda] * x[j - 1]);
            }
            rwork[i - 1] = tmp;
            anorm = max(anorm, tmp);
        }
    } else {
        for (i = 1; i <= n; i = i + 1) {
            tmp = 0.0;
            for (j = 1; j <= i; j = j + 1) {
                tmp += abs1(a[(i - 1) + (j - 1) * lda] * x[j - 1]);
            }
            for (j = i + 1; j <= n; j = j + 1) {
                tmp += abs1(a[(j - 1) + (i - 1) * lda] * x[j - 1]);
            }
            rwork[i - 1] = tmp;
            anorm = max(anorm, tmp);
        }
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return_value = 1.0;
        return return_value;
    } else if (anorm == 0.0) {
        return return_value;
    }
    //
    //     Estimate the norm of inv(op(A)).
    //
    ainvnm = 0.0;
    //
    kase = 0;
statement_10:
    Clacn2(n, &work[(n + 1) - 1], work, ainvnm, kase, isave);
    if (kase != 0) {
        if (kase == 2) {
            //
            //           Multiply by R.
            //
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = work[i - 1] * rwork[i - 1];
            }
            //
            if (up) {
                Cpotrs("U", n, 1, af, ldaf, work, n, info);
            } else {
                Cpotrs("L", n, 1, af, ldaf, work, n, info);
            }
            //
            //           Multiply by inv(X).
            //
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = work[i - 1] / x[i - 1];
            }
        } else {
            //
            //           Multiply by inv(X**H).
            //
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = work[i - 1] / x[i - 1];
            }
            //
            if (up) {
                Cpotrs("U", n, 1, af, ldaf, work, n, info);
            } else {
                Cpotrs("L", n, 1, af, ldaf, work, n, info);
            }
            //
            //           Multiply by R.
            //
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = work[i - 1] * rwork[i - 1];
            }
        }
        goto statement_10;
    }
    //
    //     Compute the estimate of the reciprocal condition number.
    //
    if (ainvnm != 0.0) {
        return_value = 1.0 / ainvnm;
    }
    //
    return return_value;
    //
}
