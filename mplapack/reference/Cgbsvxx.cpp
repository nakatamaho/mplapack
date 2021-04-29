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

void Cgbsvxx(const char *fact, const char *trans, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, COMPLEX *ab, INTEGER const ldab, COMPLEX *afb, INTEGER const ldafb, INTEGER *ipiv, char *equed, REAL *r, REAL *c, COMPLEX *b, INTEGER const ldb, COMPLEX *x, INTEGER const ldx, REAL rcond, REAL &rpvgrw, REAL *berr, INTEGER const n_err_bnds, REAL *err_bnds_norm, REAL *err_bnds_comp, INTEGER const nparams, REAL *params, COMPLEX *work, REAL *rwork, INTEGER &info) {
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
    bool notran = Mlsame(trans, "N");
    REAL smlnum = Rlamch("Safe minimum");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    bool rowequ = false;
    bool colequ = false;
    if (nofact || equil) {
        *equed = 'N';
        rowequ = false;
        colequ = false;
    } else {
        rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
        colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
    }
    //
    //     Default is failure.  If an input parameter is wrong or
    //     factorization fails, make everything look horrible.  Only the
    //     pivot growth is set here, the rest is initialized in Cgbrfsx.
    //
    const REAL zero = 0.0;
    rpvgrw = zero;
    //
    //     Test the input parameters.  PARAMS is not tested until RgerFSX.
    //
    REAL rcmin = 0.0;
    REAL rcmax = 0.0;
    INTEGER j = 0;
    REAL rowcnd = 0.0;
    REAL colcnd = 0.0;
    if (!nofact && !equil && !Mlsame(fact, "F")) {
        info = -1;
    } else if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (kl < 0) {
        info = -4;
    } else if (ku < 0) {
        info = -5;
    } else if (nrhs < 0) {
        info = -6;
    } else if (ldab < kl + ku + 1) {
        info = -8;
    } else if (ldafb < 2 * kl + ku + 1) {
        info = -10;
    } else if (Mlsame(fact, "F") && !(rowequ || colequ || Mlsame(equed, "N"))) {
        info = -12;
    } else {
        if (rowequ) {
            rcmin = bignum;
            rcmax = zero;
            for (j = 1; j <= n; j = j + 1) {
                rcmin = min(rcmin, r[j - 1]);
                rcmax = max(rcmax, r[j - 1]);
            }
            if (rcmin <= zero) {
                info = -13;
            } else if (n > 0) {
                rowcnd = max(rcmin, smlnum) / min(rcmax, bignum);
            } else {
                rowcnd = one;
            }
        }
        if (colequ && info == 0) {
            rcmin = bignum;
            rcmax = zero;
            for (j = 1; j <= n; j = j + 1) {
                rcmin = min(rcmin, c[j - 1]);
                rcmax = max(rcmax, c[j - 1]);
            }
            if (rcmin <= zero) {
                info = -14;
            } else if (n > 0) {
                colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
            } else {
                colcnd = one;
            }
        }
        if (info == 0) {
            if (ldb < max((INTEGER)1, n)) {
                info = -15;
            } else if (ldx < max((INTEGER)1, n)) {
                info = -16;
            }
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgbsvxx", -info);
        return;
    }
    //
    REAL amax = 0.0;
    INTEGER infequ = 0;
    if (equil) {
        //
        //     Compute row and column scalings to equilibrate the matrix A.
        //
        Cgbequb(n, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, infequ);
        if (infequ == 0) {
            //
            //     Equilibrate the matrix.
            //
            Claqgb(n, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
            rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
            colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
        }
        //
        //     If the scaling factors are not applied, set them to 1.0.
        //
        if (!rowequ) {
            for (j = 1; j <= n; j = j + 1) {
                r[j - 1] = 1.0;
            }
        }
        if (!colequ) {
            for (j = 1; j <= n; j = j + 1) {
                c[j - 1] = 1.0;
            }
        }
    }
    //
    //     Scale the right-hand side.
    //
    if (notran) {
        if (rowequ) {
            Clascl2(n, nrhs, r, b, ldb);
        }
    } else {
        if (colequ) {
            Clascl2(n, nrhs, c, b, ldb);
        }
    }
    //
    INTEGER i = 0;
    if (nofact || equil) {
        //
        //        Compute the LU factorization of A.
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = kl + 1; i <= 2 * kl + ku + 1; i = i + 1) {
                afb[(i - 1) + (j - 1) * ldafb] = ab[((i - kl) - 1) + (j - 1) * ldab];
            }
        }
        Cgbtrf(n, n, kl, ku, afb, ldafb, ipiv, info);
        //
        //        Return if INFO is non-zero.
        //
        if (info > 0) {
            //
            //           Pivot in column INFO is exactly 0
            //           Compute the reciprocal pivot growth factor of the
            //           leading rank-deficient INFO columns of A.
            //
            rpvgrw = Cla_gbrpvgrw(n, kl, ku, info, ab, ldab, afb, ldafb);
            return;
        }
    }
    //
    //     Compute the reciprocal pivot growth factor RPVGRW.
    //
    rpvgrw = Cla_gbrpvgrw(n, kl, ku, n, ab, ldab, afb, ldafb);
    //
    //     Compute the solution matrix X.
    //
    Clacpy("Full", n, nrhs, b, ldb, x, ldx);
    Cgbtrs(trans, n, kl, ku, nrhs, afb, ldafb, ipiv, x, ldx, info);
    //
    //     Use iterative refinement to improve the computed solution and
    //     compute error bounds and backward error estimates for it.
    //
    Cgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    //
    //     Scale solutions.
    //
    if (colequ && notran) {
        Clascl2(n, nrhs, c, x, ldx);
    } else if (rowequ && !notran) {
        Clascl2(n, nrhs, r, x, ldx);
    }
    //
    //     End of Cgbsvxx
    //
}
