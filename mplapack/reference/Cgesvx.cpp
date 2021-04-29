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

void Cgesvx(const char *fact, const char *trans, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *af, INTEGER const ldaf, INTEGER *ipiv, char *equed, REAL *r, REAL *c, COMPLEX *b, INTEGER const ldb, COMPLEX *x, INTEGER const ldx, REAL &rcond, REAL *ferr, REAL *berr, COMPLEX *work, REAL *rwork, INTEGER &info) {
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
    info = 0;
    bool nofact = Mlsame(fact, "N");
    bool equil = Mlsame(fact, "E");
    bool notran = Mlsame(trans, "N");
    bool rowequ = false;
    bool colequ = false;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    if (nofact || equil) {
        *equed = 'N';
        rowequ = false;
        colequ = false;
    } else {
        rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
        colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
        smlnum = Rlamch("Safe minimum");
        bignum = one / smlnum;
    }
    //
    //     Test the input parameters.
    //
    REAL rcmin = 0.0;
    const REAL zero = 0.0;
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
    } else if (nrhs < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldaf < max((INTEGER)1, n)) {
        info = -8;
    } else if (Mlsame(fact, "F") && !(rowequ || colequ || Mlsame(equed, "N"))) {
        info = -10;
    } else {
        if (rowequ) {
            rcmin = bignum;
            rcmax = zero;
            for (j = 1; j <= n; j = j + 1) {
                rcmin = min(rcmin, r[j - 1]);
                rcmax = max(rcmax, r[j - 1]);
            }
            if (rcmin <= zero) {
                info = -11;
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
                info = -12;
            } else if (n > 0) {
                colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
            } else {
                colcnd = one;
            }
        }
        if (info == 0) {
            if (ldb < max((INTEGER)1, n)) {
                info = -14;
            } else if (ldx < max((INTEGER)1, n)) {
                info = -16;
            }
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgesvx", -info);
        return;
    }
    //
    REAL amax = 0.0;
    INTEGER infequ = 0;
    if (equil) {
        //
        //        Compute row and column scalings to equilibrate the matrix A.
        //
        Cgeequ(n, n, a, lda, r, c, rowcnd, colcnd, amax, infequ);
        if (infequ == 0) {
            //
            //           Equilibrate the matrix.
            //
            Claqge(n, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
            rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
            colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
        }
    }
    //
    //     Scale the right hand side.
    //
    INTEGER i = 0;
    if (notran) {
        if (rowequ) {
            for (j = 1; j <= nrhs; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = r[i - 1] * b[(i - 1) + (j - 1) * ldb];
                }
            }
        }
    } else if (colequ) {
        for (j = 1; j <= nrhs; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = c[i - 1] * b[(i - 1) + (j - 1) * ldb];
            }
        }
    }
    //
    REAL rpvgrw = 0.0;
    if (nofact || equil) {
        //
        //        Compute the LU factorization of A.
        //
        Clacpy("Full", n, n, a, lda, af, ldaf);
        Cgetrf(n, n, af, ldaf, ipiv, info);
        //
        //        Return if INFO is non-zero.
        //
        if (info > 0) {
            //
            //           Compute the reciprocal pivot growth factor of the
            //           leading rank-deficient INFO columns of A.
            //
            rpvgrw = Clantr("M", "U", "N", info, info, af, ldaf, rwork);
            if (rpvgrw == zero) {
                rpvgrw = one;
            } else {
                rpvgrw = Clange("M", n, info, a, lda, rwork) / rpvgrw;
            }
            rwork[1 - 1] = rpvgrw;
            rcond = zero;
            return;
        }
    }
    //
    //     Compute the norm of the matrix A and the
    //     reciprocal pivot growth factor RPVGRW.
    //
    char norm;
    if (notran) {
        norm = "1";
    } else {
        norm = "I";
    }
    REAL anorm = Clange(norm, n, n, a, lda, rwork);
    rpvgrw = Clantr("M", "U", "N", n, n, af, ldaf, rwork);
    if (rpvgrw == zero) {
        rpvgrw = one;
    } else {
        rpvgrw = Clange("M", n, n, a, lda, rwork) / rpvgrw;
    }
    //
    //     Compute the reciprocal of the condition number of A.
    //
    Cgecon(norm, n, af, ldaf, anorm, rcond, work, rwork, info);
    //
    //     Compute the solution matrix X.
    //
    Clacpy("Full", n, nrhs, b, ldb, x, ldx);
    Cgetrs(trans, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    //
    //     Use iterative refinement to improve the computed solution and
    //     compute error bounds and backward error estimates for it.
    //
    Cgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    //
    //     Transform the solution matrix X to a solution of the original
    //     system.
    //
    if (notran) {
        if (colequ) {
            for (j = 1; j <= nrhs; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    x[(i - 1) + (j - 1) * ldx] = c[i - 1] * x[(i - 1) + (j - 1) * ldx];
                }
            }
            for (j = 1; j <= nrhs; j = j + 1) {
                ferr[j - 1] = ferr[j - 1] / colcnd;
            }
        }
    } else if (rowequ) {
        for (j = 1; j <= nrhs; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                x[(i - 1) + (j - 1) * ldx] = r[i - 1] * x[(i - 1) + (j - 1) * ldx];
            }
        }
        for (j = 1; j <= nrhs; j = j + 1) {
            ferr[j - 1] = ferr[j - 1] / rowcnd;
        }
    }
    //
    //     Set INFO = N+1 if the matrix is singular to working precision.
    //
    if (rcond < Rlamch("Epsilon")) {
        info = n + 1;
    }
    //
    rwork[1 - 1] = rpvgrw;
    //
    //     End of Cgesvx
    //
}
