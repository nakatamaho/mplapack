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

void Rsposv(const char *uplo, INTEGER const &n, INTEGER const &nrhs, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *x, INTEGER const &ldx, REAL *work, arr_cref<float> swork, INTEGER &iter, INTEGER &info) {
    const bool doitref = true;
    REAL anrm = 0.0;
    REAL eps = 0.0;
    const REAL bwdmax = 1.00f;
    REAL cte = 0.0;
    INTEGER ptsa = 0;
    INTEGER ptsx = 0;
    const REAL negone = -1.0;
    const REAL one = 1.0;
    INTEGER i = 0;
    REAL xnrm = 0.0;
    REAL rnrm = 0.0;
    INTEGER iiter = 0;
    const INTEGER itermax = 30;
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
    //
    //     .. Local Scalars ..
    //
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    iter = 0;
    //
    //     Test the input parameters.
    //
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldx < max((INTEGER)1, n)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rsposv", -info);
        return;
    }
    //
    //     Quick return if (N.EQ.0).
    //
    if (n == 0) {
        return;
    }
    //
    //     Skip single precision iterative refinement if a priori slower
    //     than REAL precision factorization.
    //
    if (!doitref) {
        iter = -1;
        goto statement_40;
    }
    //
    //     Compute some constants.
    //
    anrm = Rlansy[("I" - 1) + (uplo - 1) * ldRlansy];
    eps = dlamch("Epsilon");
    cte = anrm * eps * sqrt(n.real()) * bwdmax;
    //
    //     Set the indices PTSA, PTSX for referencing SA and SX in SWORK.
    //
    ptsa = 1;
    ptsx = ptsa + n * n;
    //
    //     Convert B from REAL precision to single precision and store the
    //     result in SX.
    //
    Rlag2s(n, nrhs, b, ldb, swork[ptsx - 1], n, info);
    //
    if (info != 0) {
        iter = -2;
        goto statement_40;
    }
    //
    //     Convert A from REAL precision to single precision and store the
    //     result in SA.
    //
    Rlat2s(uplo, n, a, lda, swork[ptsa - 1], n, info);
    //
    if (info != 0) {
        iter = -2;
        goto statement_40;
    }
    //
    //     Compute the Cholesky factorization of SA.
    //
    spotrf(uplo, n, swork[ptsa - 1], n, info);
    //
    if (info != 0) {
        iter = -3;
        goto statement_40;
    }
    //
    //     Solve the system SA*SX = SB.
    //
    spotrs(uplo, n, nrhs, swork[ptsa - 1], n, swork[ptsx - 1], n, info);
    //
    //     Convert SX back to REAL precision
    //
    slag2d(n, nrhs, swork[ptsx - 1], n, x, ldx, info);
    //
    //     Compute R = B - AX (R is WORK).
    //
    Rlacpy("All", n, nrhs, b, ldb, work, n);
    //
    Rsymm("Left", uplo, n, nrhs, negone, a, lda, x, ldx, one, work, n);
    //
    //     Check whether the NRHS normwise backward errors satisfy the
    //     stopping criterion. If yes, set ITER=0 and return.
    //
    for (i = 1; i <= nrhs; i = i + 1) {
        xnrm = abs(x[(iRamax[(n - 1) + (x[(i - 1) * ldx] - 1) * ldiRamax] - 1) + (i - 1) * ldx]);
        rnrm = abs(work[(iRamax[(n - 1) + (work[(i - 1) * ldwork] - 1) * ldiRamax] - 1) + (i - 1) * ldwork]);
        if (rnrm > xnrm * cte) {
            goto statement_10;
        }
    }
    //
    //     If we are here, the NRHS normwise backward errors satisfy the
    //     stopping criterion. We are good to exit.
    //
    iter = 0;
    return;
//
statement_10:
    //
    for (iiter = 1; iiter <= itermax; iiter = iiter + 1) {
        //
        //        Convert R (in WORK) from REAL precision to single precision
        //        and store the result in SX.
        //
        Rlag2s(n, nrhs, work, n, swork[ptsx - 1], n, info);
        //
        if (info != 0) {
            iter = -2;
            goto statement_40;
        }
        //
        //        Solve the system SA*SX = SR.
        //
        spotrs(uplo, n, nrhs, swork[ptsa - 1], n, swork[ptsx - 1], n, info);
        //
        //        Convert SX back to REAL precision and update the current
        //        iterate.
        //
        slag2d(n, nrhs, swork[ptsx - 1], n, work, n, info);
        //
        for (i = 1; i <= nrhs; i = i + 1) {
            Raxpy(n, one, work[(i - 1) * ldwork], 1, x[(i - 1) * ldx], 1);
        }
        //
        //        Compute R = B - AX (R is WORK).
        //
        Rlacpy("All", n, nrhs, b, ldb, work, n);
        //
        Rsymm("L", uplo, n, nrhs, negone, a, lda, x, ldx, one, work, n);
        //
        //        Check whether the NRHS normwise backward errors satisfy the
        //        stopping criterion. If yes, set ITER=IITER>0 and return.
        //
        for (i = 1; i <= nrhs; i = i + 1) {
            xnrm = abs(x[(iRamax[(n - 1) + (x[(i - 1) * ldx] - 1) * ldiRamax] - 1) + (i - 1) * ldx]);
            rnrm = abs(work[(iRamax[(n - 1) + (work[(i - 1) * ldwork] - 1) * ldiRamax] - 1) + (i - 1) * ldwork]);
            if (rnrm > xnrm * cte) {
                goto statement_20;
            }
        }
        //
        //        If we are here, the NRHS normwise backward errors satisfy the
        //        stopping criterion, we are good to exit.
        //
        iter = iiter;
        //
        return;
    //
    statement_20:;
        //
    }
    //
    //     If we are at this place of the code, this is because we have
    //     performed ITER=ITERMAX iterations and never satisfied the
    //     stopping criterion, set up the ITER flag accordingly and follow
    //     up on REAL precision routine.
    //
    iter = -itermax - 1;
//
statement_40:
    //
    //     Single-precision iterative refinement failed to converge to a
    //     satisfactory solution, so we resort to REAL precision.
    //
    Rpotrf(uplo, n, a, lda, info);
    //
    if (info != 0) {
        return;
    }
    //
    Rlacpy("All", n, nrhs, b, ldb, x, ldx);
    Rpotrs(uplo, n, nrhs, a, lda, x, ldx, info);
    //
    //     End of Rsposv.
    //
}
