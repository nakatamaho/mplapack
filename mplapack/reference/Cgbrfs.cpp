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

inline REAL abs1(COMPLEX zdum) {
   return abs(zdum.real()) + abs(zdum.imag());
}

void Cgbrfs(const char *trans, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, COMPLEX *ab, INTEGER const ldab, COMPLEX *afb, INTEGER const ldafb, INTEGER *ipiv, COMPLEX *b, INTEGER const ldb, COMPLEX *x, INTEGER const ldx, REAL *ferr, REAL *berr, COMPLEX *work, REAL *rwork, INTEGER &info) {
    COMPLEX zdum = 0.0;
    bool notran = false;
    INTEGER j = 0;
    const REAL zero = 0.0;
    char transn;
    char transt;
    INTEGER nz = 0;
    REAL eps = 0.0;
    REAL safmin = 0.0;
    REAL safe1 = 0.0;
    REAL safe2 = 0.0;
    INTEGER count = 0;
    const REAL three = 3.0e+0;
    REAL lstres = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER i = 0;
    INTEGER k = 0;
    INTEGER kk = 0;
    REAL xk = 0.0;
    REAL s = 0.0;
    const REAL two = 2.0e+0;
    const INTEGER itmax = 5;
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1(zdum) = abs(zdum.real()) + abs(zdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0) {
        info = -3;
    } else if (ku < 0) {
        info = -4;
    } else if (nrhs < 0) {
        info = -5;
    } else if (ldab < kl + ku + 1) {
        info = -7;
    } else if (ldafb < 2 * kl + ku + 1) {
        info = -9;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -12;
    } else if (ldx < max((INTEGER)1, n)) {
        info = -14;
    }
    if (info != 0) {
        Mxerbla("Cgbrfs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        for (j = 1; j <= nrhs; j = j + 1) {
            ferr[j - 1] = zero;
            berr[j - 1] = zero;
        }
        return;
    }
    //
    if (notran) {
        transn = 'N';
        transt = 'C';
    } else {
        transn = 'C';
        transt = 'N';
    }
    //
    //     NZ = maximum number of nonzero elements in each row of A, plus 1
    //
    nz = min(kl + ku + 2, n + 1);
    eps = Rlamch("Epsilon");
    safmin = Rlamch("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
    //
    //     Do for each right hand side
    //
    for (j = 1; j <= nrhs; j = j + 1) {
        //
        count = 1;
        lstres = three;
    statement_20:
        //
        //        Loop until stopping criterion is satisfied.
        //
        //        Compute residual R = B - op(A) * X,
        //        where op(A) = A, A**T, or A**H, depending on TRANS.
        //
        Ccopy(n, &b[(j - 1) * ldb], 1, work, 1);
        Cgbmv(trans, n, n, kl, ku, -cone, ab, ldab, &x[(j - 1) * ldx], 1, cone, work, 1);
        //
        //        Compute componentwise relative backward error from formula
        //
        //        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
        //
        //        where abs(Z) is the componentwise absolute value of the matrix
        //        or vector Z.  If the i-th component of the denominator is less
        //        than SAFE2, then SAFE1 is added to the i-th components of the
        //        numerator and denominator before dividing.
        //
        for (i = 1; i <= n; i = i + 1) {
            rwork[i - 1] = abs1(b[(i - 1) + (j - 1) * ldb]);
        }
        //
        //        Compute abs(op(A))*abs(X) + abs(B).
        //
        if (notran) {
            for (k = 1; k <= n; k = k + 1) {
                kk = ku + 1 - k;
                xk = abs1(x[(k - 1) + (j - 1) * ldx]);
                for (i = max((INTEGER)1, k - ku); i <= min(n, k + kl); i = i + 1) {
                    rwork[i - 1] += abs1(ab[((kk + i) - 1) + (k - 1) * ldab]) * xk;
                }
            }
        } else {
            for (k = 1; k <= n; k = k + 1) {
                s = zero;
                kk = ku + 1 - k;
                for (i = max((INTEGER)1, k - ku); i <= min(n, k + kl); i = i + 1) {
                    s += abs1(ab[((kk + i) - 1) + (k - 1) * ldab]) * abs1(x[(i - 1) + (j - 1) * ldx]);
                }
                rwork[k - 1] += s;
            }
        }
        s = zero;
        for (i = 1; i <= n; i = i + 1) {
            if (rwork[i - 1] > safe2) {
                s = max(s, abs1(work[i - 1]) / rwork[i - 1]);
            } else {
                s = max(s, (abs1(work[i - 1]) + safe1) / (rwork[i - 1] + safe1));
            }
        }
        berr[j - 1] = s;
        //
        //        Test stopping criterion. Continue iterating if
        //           1) The residual BERR(J) is larger than machine epsilon, and
        //           2) BERR(J) decreased by at least a factor of 2 during the
        //              last iteration, and
        //           3) At most ITMAX iterations tried.
        //
        if (berr[j - 1] > eps && two * berr[j - 1] <= lstres && count <= itmax) {
            //
            //           Update solution and try again.
            //
            Cgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
            Caxpy(n, cone, work, 1, &x[(j - 1) * ldx], 1);
            lstres = berr[j - 1];
            count++;
            goto statement_20;
        }
        //
        //        Bound error from formula
        //
        //        norm(X - XTRUE) / norm(X) .le. FERR =
        //        norm( abs(inv(op(A)))*
        //           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
        //
        //        where
        //          norm(Z) is the magnitude of the largest component of Z
        //          inv(op(A)) is the inverse of op(A)
        //          abs(Z) is the componentwise absolute value of the matrix or
        //             vector Z
        //          NZ is the maximum number of nonzeros in any row of A, plus 1
        //          EPS is machine epsilon
        //
        //        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
        //        is incremented by SAFE1 if the i-th component of
        //        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
        //
        //        Use Clacn2 to estimate the infinity-norm of the matrix
        //           inv(op(A)) * diag(W),
        //        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
        //
        for (i = 1; i <= n; i = i + 1) {
            if (rwork[i - 1] > safe2) {
                rwork[i - 1] = abs1(work[i - 1]) + nz * eps * rwork[i - 1];
            } else {
                rwork[i - 1] = abs1(work[i - 1]) + nz * eps * rwork[i - 1] + safe1;
            }
        }
        //
        kase = 0;
    statement_100:
        Clacn2(n, &work[(n + 1) - 1], work, ferr[j - 1], kase, isave);
        if (kase != 0) {
            if (kase == 1) {
                //
                //              Multiply by diag(W)*inv(op(A)**H).
                //
                Cgbtrs(&transt, n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = rwork[i - 1] * work[i - 1];
                }
            } else {
                //
                //              Multiply by inv(op(A))*diag(W).
                //
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = rwork[i - 1] * work[i - 1];
                }
                Cgbtrs(&transn, n, kl, ku, 1, afb, ldafb, ipiv, work, n, info);
            }
            goto statement_100;
        }
        //
        //        Normalize error.
        //
        lstres = zero;
        for (i = 1; i <= n; i = i + 1) {
            lstres = max(lstres, abs1(x[(i - 1) + (j - 1) * ldx]));
        }
        if (lstres != zero) {
            ferr[j - 1] = ferr[j - 1] / lstres;
        }
        //
    }
    //
    //     End of Cgbrfs
    //
}
