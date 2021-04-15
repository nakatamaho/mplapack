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

void Rgtcon(const char *norm, INTEGER const n, REAL *dl, REAL *d, REAL *du, REAL *du2, INTEGER *ipiv, REAL const anorm, REAL &rcond, REAL *work, INTEGER *iwork, INTEGER &info) {
    bool onenrm = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER i = 0;
    REAL ainvnm = 0.0;
    INTEGER kase1 = 0;
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
    //     .. Executable Statements ..
    //
    //     Test the input arguments.
    //
    info = 0;
    onenrm = norm == "1" || Mlsame(norm, "O");
    if (!onenrm && !Mlsame(norm, "I")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (anorm < zero) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Rgtcon", -info);
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
    //     Check that D(1:N) is non-zero.
    //
    for (i = 1; i <= n; i = i + 1) {
        if (d[i - 1] == zero) {
            return;
        }
    }
    //
    ainvnm = zero;
    if (onenrm) {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kase = 0;
statement_20:
    Rlacn2(n, &work[(n + 1) - 1], work, iwork, ainvnm, kase, isave);
    if (kase != 0) {
        if (kase == kase1) {
            //
            //           Multiply by inv(U)*inv(L).
            //
            Rgttrs("No transpose", n, 1, dl, d, du, du2, ipiv, work, n, info);
        } else {
            //
            //           Multiply by inv(L**T)*inv(U**T).
            //
            Rgttrs("Transpose", n, 1, dl, d, du, du2, ipiv, work, n, info);
        }
        goto statement_20;
    }
    //
    //     Compute the estimate of the reciprocal condition number.
    //
    if (ainvnm != zero) {
        rcond = (one / ainvnm) / anorm;
    }
    //
    //     End of Rgtcon
    //
}
