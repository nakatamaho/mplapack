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

void Cposvxx(const char *fact, const char *uplo, INTEGER const &n, INTEGER const &nrhs, COMPLEX *a, INTEGER const &lda, COMPLEX *af, INTEGER const &ldaf, str_ref equed, REAL *s, COMPLEX *b, INTEGER const &ldb, COMPLEX *x, INTEGER const &ldx, REAL const &rcond, REAL &rpvgrw, REAL *berr, INTEGER const &n_err_bnds, REAL *err_bnds_norm, REAL *err_bnds_comp, INTEGER const &nparams, REAL *params, COMPLEX *work, REAL *rwork, INTEGER &info) {
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
    //  ==================================================================
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
    info = 0;
    bool nofact = Mlsame(fact, "N");
    bool equil = Mlsame(fact, "E");
    REAL smlnum = dlamch("Safe minimum");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    bool rcequ = false;
    if (nofact || equil) {
        equed = "N";
        rcequ = false;
    } else {
        rcequ = Mlsame(equed, "Y");
    }
    //
    //     Default is failure.  If an input parameter is wrong or
    //     factorization fails, make everything look horrible.  Only the
    //     pivot growth is set here, the rest is initialized in Cporfsx.
    //
    const REAL zero = 0.0;
    rpvgrw = zero;
    //
    //     Test the input parameters.  PARAMS is not tested until Cporfsx.
    //
    REAL smin = 0.0;
    REAL smax = 0.0;
    INTEGER j = 0;
    REAL scond = 0.0;
    if (!nofact && !equil && !Mlsame(fact, "F")) {
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
    } else if (Mlsame(fact, "F") && !(rcequ || Mlsame(equed, "N"))) {
        info = -9;
    } else {
        if (rcequ) {
            smin = bignum;
            smax = zero;
            for (j = 1; j <= n; j = j + 1) {
                smin = min(smin, s[j - 1]);
                smax = max(smax, s[j - 1]);
            }
            if (smin <= zero) {
                info = -10;
            } else if (n > 0) {
                scond = max(smin, smlnum) / min(smax, bignum);
            } else {
                scond = one;
            }
        }
        if (info == 0) {
            if (ldb < max((INTEGER)1, n)) {
                info = -12;
            } else if (ldx < max((INTEGER)1, n)) {
                info = -14;
            }
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cposvxx", -info);
        return;
    }
    //
    REAL amax = 0.0;
    INTEGER infequ = 0;
    if (equil) {
        //
        //     Compute row and column scalings to equilibrate the matrix A.
        //
        Cpoequb(n, a, lda, s, scond, amax, infequ);
        if (infequ == 0) {
            //
            //     Equilibrate the matrix.
            //
            Claqhe(uplo, n, a, lda, s, scond, amax, equed);
            rcequ = Mlsame(equed, "Y");
        }
    }
    //
    //     Scale the right-hand side.
    //
    if (rcequ) {
        Clascl2(n, nrhs, s, b, ldb);
    }
    //
    if (nofact || equil) {
        //
        //        Compute the Cholesky factorization of A.
        //
        Clacpy(uplo, n, n, a, lda, af, ldaf);
        Cpotrf(uplo, n, af, ldaf, info);
        //
        //        Return if INFO is non-zero.
        //
        if (info > 0) {
            //
            //           Pivot in column INFO is exactly 0
            //           Compute the reciprocal pivot growth factor of the
            //           leading rank-deficient INFO columns of A.
            //
            rpvgrw = Cla_porpvgrw[(uplo - 1) + (n - 1) * ldCla_porpvgrw];
            return;
        }
    }
    //
    //     Compute the reciprocal pivot growth factor RPVGRW.
    //
    rpvgrw = Cla_porpvgrw[(uplo - 1) + (n - 1) * ldCla_porpvgrw];
    //
    //     Compute the solution matrix X.
    //
    Clacpy("Full", n, nrhs, b, ldb, x, ldx);
    Cpotrs(uplo, n, nrhs, af, ldaf, x, ldx, info);
    //
    //     Use iterative refinement to improve the computed solution and
    //     compute error bounds and backward error estimates for it.
    //
    Cporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    //
    //     Scale solutions.
    //
    if (rcequ) {
        Clascl2(n, nrhs, s, x, ldx);
    }
    //
    //     End of Cposvxx
    //
}
