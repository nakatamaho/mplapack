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

void Rsysvx(const char *fact, const char *uplo, INTEGER const &n, INTEGER const &nrhs, REAL *a, INTEGER const &lda, REAL *af, INTEGER const &ldaf, INTEGER *ipiv, REAL *b, INTEGER const &ldb, REAL *x, INTEGER const &ldx, REAL &rcond, REAL *ferr, REAL *berr, REAL *work, INTEGER const &lwork, INTEGER *iwork, INTEGER &info) {
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
    bool lquery = (lwork == -1);
    if (!nofact && !Mlsame(fact, "F")) {
        info = -1;
    } else if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (nrhs < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldaf < max((INTEGER)1, n)) {
        info = -8;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -11;
    } else if (ldx < max((INTEGER)1, n)) {
        info = -13;
    } else if (lwork < max((INTEGER)1, 3 * n) && !lquery) {
        info = -18;
    }
    //
    INTEGER lwkopt = 0;
    INTEGER nb = 0;
    if (info == 0) {
        lwkopt = max((INTEGER)1, 3 * n);
        if (nofact) {
            nb = iMlaenv[("Rsytrf" - 1) * ldiMlaenv];
            lwkopt = max(lwkopt, n * nb);
        }
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rsysvx", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    const REAL zero = 0.0;
    if (nofact) {
        //
        //        Compute the factorization A = U*D*U**T or A = L*D*L**T.
        //
        Rlacpy(uplo, n, n, a, lda, af, ldaf);
        Rsytrf(uplo, n, af, ldaf, ipiv, work, lwork, info);
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
    REAL anorm = Rlansy[("I" - 1) + (uplo - 1) * ldRlansy];
    //
    //     Compute the reciprocal of the condition number of A.
    //
    Rsycon(uplo, n, af, ldaf, ipiv, anorm, rcond, work, iwork, info);
    //
    //     Compute the solution vectors X.
    //
    Rlacpy("Full", n, nrhs, b, ldb, x, ldx);
    Rsytrs(uplo, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    //
    //     Use iterative refinement to improve the computed solutions and
    //     compute error bounds and backward error estimates for them.
    //
    Rsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    //
    //     Set INFO = N+1 if the matrix is singular to working precision.
    //
    if (rcond < dlamch("Epsilon")) {
        info = n + 1;
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Rsysvx
    //
}
