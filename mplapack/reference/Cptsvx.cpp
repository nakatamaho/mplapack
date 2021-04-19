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

void Cptsvx(const char *fact, INTEGER const n, INTEGER const nrhs, REAL *d, COMPLEX *e, REAL *df, COMPLEX *ef, COMPLEX *b, INTEGER const ldb, COMPLEX *x, INTEGER const ldx, REAL &rcond, REAL *ferr, REAL *berr, COMPLEX *work, REAL *rwork, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
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
    //     Test the input parameters.
    //
    info = 0;
    bool nofact = Mlsame(fact, "N");
    if (!nofact && !Mlsame(fact, "F")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if (ldx < max((INTEGER)1, n)) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Cptsvx", -info);
        return;
    }
    //
    const REAL zero = 0.0;
    if (nofact) {
        //
        //        Compute the L*D*L**H (or U**H*D*U) factorization of A.
        //
        Rcopy(n, d, 1, df, 1);
        if (n > 1) {
            Ccopy(n - 1, e, 1, ef, 1);
        }
        Cpttrf(n, df, ef, info);
        //
        //        Return if INFO is non-zero.
        //
        if (info > 0) {
            rcond = zero;
            return;
        }
    }
    //
    //     Compute the norm of the matrix A.
    //
    REAL anorm = Clanht("1", n, d, e);
    //
    //     Compute the reciprocal of the condition number of A.
    //
    Cptcon(n, df, ef, anorm, rcond, rwork, info);
    //
    //     Compute the solution vectors X.
    //
    Clacpy("Full", n, nrhs, b, ldb, x, ldx);
    Cpttrs("Lower", n, nrhs, df, ef, x, ldx, info);
    //
    //     Use iterative refinement to improve the computed solutions and
    //     compute error bounds and backward error estimates for them.
    //
    Cptrfs("Lower", n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    //
    //     Set INFO = N+1 if the matrix is singular to working precision.
    //
    if (rcond < Rlamch("Epsilon")) {
        info = n + 1;
    }
    //
    //     End of Cptsvx
    //
}
