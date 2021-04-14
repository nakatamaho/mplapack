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

void Rlasd2(INTEGER const nl, INTEGER const nr, INTEGER const sqre, INTEGER &k, REAL *d, REAL *z, REAL const alpha, REAL const beta, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, REAL *dsigma, REAL *u2, INTEGER const ldu2, REAL *vt2, INTEGER const ldvt2, INTEGER *idxp, INTEGER *idx, INTEGER *idxc, INTEGER *idxq, INTEGER *coltyp, INTEGER &info) {
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER nlp1 = 0;
    INTEGER nlp2 = 0;
    REAL z1 = 0.0;
    INTEGER i = 0;
    INTEGER idxi = 0;
    REAL eps = 0.0;
    REAL tol = 0.0;
    const REAL eight = 8.0e+0;
    INTEGER k2 = 0;
    INTEGER j = 0;
    INTEGER jprev = 0;
    REAL s = 0.0;
    REAL c = 0.0;
    REAL tau = 0.0;
    const REAL zero = 0.0;
    INTEGER idxjp = 0;
    INTEGER idxj = 0;
    INTEGER ctot[4];
    INTEGER ct = 0;
    INTEGER psm[4];
    INTEGER jp = 0;
    const REAL two = 2.0e+0;
    REAL hlftol = 0.0;
    const REAL one = 1.0;
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. Local Arrays ..
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
    //
    if (nl < 1) {
        info = -1;
    } else if (nr < 1) {
        info = -2;
    } else if ((sqre != 1) && (sqre != 0)) {
        info = -3;
    }
    //
    n = nl + nr + 1;
    m = n + sqre;
    //
    if (ldu < n) {
        info = -10;
    } else if (ldvt < m) {
        info = -12;
    } else if (ldu2 < n) {
        info = -15;
    } else if (ldvt2 < m) {
        info = -17;
    }
    if (info != 0) {
        Mxerbla("Rlasd2", -info);
        return;
    }
    //
    nlp1 = nl + 1;
    nlp2 = nl + 2;
    //
    //     Generate the first part of the vector Z; and move the singular
    //     values in the first part of D one position backward.
    //
    z1 = alpha * vt[(nlp1 - 1) + (nlp1 - 1) * ldvt];
    z[1 - 1] = z1;
    for (i = nl; i >= 1; i = i - 1) {
        z[(i + 1) - 1] = alpha * vt[(i - 1) + (nlp1 - 1) * ldvt];
        d[(i + 1) - 1] = d[i - 1];
        idxq[(i + 1) - 1] = idxq[i - 1] + 1;
    }
    //
    //     Generate the second part of the vector Z.
    //
    for (i = nlp2; i <= m; i = i + 1) {
        z[i - 1] = beta * vt[(i - 1) + (nlp2 - 1) * ldvt];
    }
    //
    //     Initialize some reference arrays.
    //
    for (i = 2; i <= nlp1; i = i + 1) {
        coltyp[i - 1] = 1;
    }
    for (i = nlp2; i <= n; i = i + 1) {
        coltyp[i - 1] = 2;
    }
    //
    //     Sort the singular values into increasing order
    //
    for (i = nlp2; i <= n; i = i + 1) {
        idxq[i - 1] += nlp1;
    }
    //
    //     DSIGMA, IDXC, IDXC, and the first column of U2
    //     are used as storage space.
    //
    for (i = 2; i <= n; i = i + 1) {
        dsigma[i - 1] = d[idxq[i - 1] - 1];
        u2[(i - 1)] = z[idxq[i - 1] - 1];
        idxc[i - 1] = coltyp[idxq[i - 1] - 1];
    }
    //
    Rlamrg(nl, nr, &dsigma[2 - 1], 1, 1, &idx[2 - 1]);
    //
    for (i = 2; i <= n; i = i + 1) {
        idxi = 1 + idx[i - 1];
        d[i - 1] = dsigma[idxi - 1];
        z[i - 1] = u2[(idxi - 1)];
        coltyp[i - 1] = idxc[idxi - 1];
    }
    //
    //     Calculate the allowable deflation tolerance
    //
    eps = Rlamch("Epsilon");
    tol = max(abs(alpha), abs(beta));
    tol = eight * eps * max(abs(d[n - 1]), tol);
    //
    //     There are 2 kinds of deflation -- first a value in the z-vector
    //     is small, second two (or more) singular values are very close
    //     together (their difference is small).
    //
    //     If the value in the z-vector is small, we simply permute the
    //     array so that the corresponding singular value is moved to the
    //     end.
    //
    //     If two values in the D-vector are close, we perform a two-sided
    //     rotation designed to make one of the corresponding z-vector
    //     entries zero, and then permute the array so that the deflated
    //     singular value is moved to the end.
    //
    //     If there are multiple singular values then the problem deflates.
    //     Here the number of equal singular values are found.  As each equal
    //     singular value is found, an elementary reflector is computed to
    //     rotate the corresponding singular subspace so that the
    //     corresponding components of Z are zero in this new basis.
    //
    k = 1;
    k2 = n + 1;
    for (j = 2; j <= n; j = j + 1) {
        if (abs(z[j - 1]) <= tol) {
            //
            //           Deflate due to small z component.
            //
            k2 = k2 - 1;
            idxp[k2 - 1] = j;
            coltyp[j - 1] = 4;
            if (j == n) {
                goto statement_120;
            }
        } else {
            jprev = j;
            goto statement_90;
        }
    }
statement_90:
    j = jprev;
statement_100:
    j++;
    if (j > n) {
        goto statement_110;
    }
    if (abs(z[j - 1]) <= tol) {
        //
        //        Deflate due to small z component.
        //
        k2 = k2 - 1;
        idxp[k2 - 1] = j;
        coltyp[j - 1] = 4;
    } else {
        //
        //        Check if singular values are close enough to allow deflation.
        //
        if (abs(d[j - 1] - d[jprev - 1]) <= tol) {
            //
            //           Deflation is possible.
            //
            s = z[jprev - 1];
            c = z[j - 1];
            //
            //           Find sqrt(a**2+b**2) without overflow or
            //           destructive underflow.
            //
            tau = Rlapy2(c, s);
            c = c / tau;
            s = -s / tau;
            z[j - 1] = tau;
            z[jprev - 1] = zero;
            //
            //           Apply back the Givens rotation to the left and right
            //           singular vector matrices.
            //
            idxjp = idxq[(idx[jprev - 1] + 1) - 1];
            idxj = idxq[(idx[j - 1] + 1) - 1];
            if (idxjp <= nlp1) {
                idxjp = idxjp - 1;
            }
            if (idxj <= nlp1) {
                idxj = idxj - 1;
            }
            Rrot(n, &u[(idxjp - 1) * ldu], 1, &u[(idxj - 1) * ldu], 1, c, s);
            Rrot(m, &vt[(idxjp - 1)], ldvt, &vt[(idxj - 1)], ldvt, c, s);
            if (coltyp[j - 1] != coltyp[jprev - 1]) {
                coltyp[j - 1] = 3;
            }
            coltyp[jprev - 1] = 4;
            k2 = k2 - 1;
            idxp[k2 - 1] = jprev;
            jprev = j;
        } else {
            k++;
            u2[(k - 1)] = z[jprev - 1];
            dsigma[k - 1] = d[jprev - 1];
            idxp[k - 1] = jprev;
            jprev = j;
        }
    }
    goto statement_100;
statement_110:
    //
    //     Record the last singular value.
    //
    k++;
    u2[(k - 1)] = z[jprev - 1];
    dsigma[k - 1] = d[jprev - 1];
    idxp[k - 1] = jprev;
//
statement_120:
    //
    //     Count up the total number of the various types of columns, then
    //     form a permutation which positions the four column types into
    //     four groups of uniform structure (although one or more of these
    //     groups may be empty).
    //
    for (j = 1; j <= 4; j = j + 1) {
        ctot[j - 1] = 0;
    }
    for (j = 2; j <= n; j = j + 1) {
        ct = coltyp[j - 1];
        ctot[ct - 1]++;
    }
    //
    //     PSM(*) = Position in SubMatrix (of types 1 through 4)
    //
    psm[1 - 1] = 2;
    psm[2 - 1] = 2 + ctot[1 - 1];
    psm[3 - 1] = psm[2 - 1] + ctot[2 - 1];
    psm[4 - 1] = psm[3 - 1] + ctot[3 - 1];
    //
    //     Fill out the IDXC array so that the permutation which it induces
    //     will place all type-1 columns first, all type-2 columns next,
    //     then all type-3's, and finally all type-4's, starting from the
    //     second column. This applies similarly to the rows of VT.
    //
    for (j = 2; j <= n; j = j + 1) {
        jp = idxp[j - 1];
        ct = coltyp[jp - 1];
        idxc[psm[ct - 1] - 1] = j;
        psm[ct - 1]++;
    }
    //
    //     Sort the singular values and corresponding singular vectors into
    //     DSIGMA, U2, and VT2 respectively.  The singular values/vectors
    //     which were not deflated go into the first K slots of DSIGMA, U2,
    //     and VT2 respectively, while those which were deflated go into the
    //     last N - K slots, except that the first column/row will be treated
    //     separately.
    //
    for (j = 2; j <= n; j = j + 1) {
        jp = idxp[j - 1];
        dsigma[j - 1] = d[jp - 1];
        idxj = idxq[(idx[idxp[idxc[j - 1] - 1] - 1] + 1) - 1];
        if (idxj <= nlp1) {
            idxj = idxj - 1;
        }
        Rcopy(n, &u[(idxj - 1) * ldu], 1, &u2[(j - 1) * ldu2], 1);
        Rcopy(m, &vt[(idxj - 1)], ldvt, &vt2[(j - 1)], ldvt2);
    }
    //
    //     Determine DSIGMA(1), DSIGMA(2) and Z(1)
    //
    dsigma[1 - 1] = zero;
    hlftol = tol / two;
    if (abs(dsigma[2 - 1]) <= hlftol) {
        dsigma[2 - 1] = hlftol;
    }
    if (m > n) {
        z[1 - 1] = Rlapy2(z1, z[m - 1]);
        if (z[1 - 1] <= tol) {
            c = one;
            s = zero;
            z[1 - 1] = tol;
        } else {
            c = z1 / z[1 - 1];
            s = z[m - 1] / z[1 - 1];
        }
    } else {
        if (abs(z1) <= tol) {
            z[1 - 1] = tol;
        } else {
            z[1 - 1] = z1;
        }
    }
    //
    //     Move the rest of the updating row to Z.
    //
    Rcopy(k - 1, &u2[(2 - 1)], 1, &z[2 - 1], 1);
    //
    //     Determine the first column of U2, the first row of VT2 and the
    //     last row of VT.
    //
    Rlaset("A", n, 1, zero, zero, u2, ldu2);
    u2[(nlp1 - 1)] = one;
    if (m > n) {
        for (i = 1; i <= nlp1; i = i + 1) {
            vt[(m - 1) + (i - 1) * ldvt] = -s * vt[(nlp1 - 1) + (i - 1) * ldvt];
            vt2[(i - 1) * ldvt2] = c * vt[(nlp1 - 1) + (i - 1) * ldvt];
        }
        for (i = nlp2; i <= m; i = i + 1) {
            vt2[(i - 1) * ldvt2] = s * vt[(m - 1) + (i - 1) * ldvt];
            vt[(m - 1) + (i - 1) * ldvt] = c * vt[(m - 1) + (i - 1) * ldvt];
        }
    } else {
        Rcopy(m, &vt[(nlp1 - 1)], ldvt, &vt2[(1 - 1)], ldvt2);
    }
    if (m > n) {
        Rcopy(m, &vt[(m - 1)], ldvt, &vt2[(m - 1)], ldvt2);
    }
    //
    //     The deflated singular values and their corresponding vectors go
    //     into the back of D, U, and V respectively.
    //
    if (n > k) {
        Rcopy(n - k, &dsigma[(k + 1) - 1], 1, &d[(k + 1) - 1], 1);
        Rlacpy("A", n, n - k, &u2[((k + 1) - 1) * ldu2], ldu2, &u[((k + 1) - 1) * ldu], ldu);
        Rlacpy("A", n - k, m, &vt2[((k + 1) - 1)], ldvt2, &vt[((k + 1) - 1)], ldvt);
    }
    //
    //     Copy CTOT into COLTYP for referencing in Rlasd3.
    //
    for (j = 1; j <= 4; j = j + 1) {
        coltyp[j - 1] = ctot[j - 1];
    }
    //
    //     End of Rlasd2
    //
}
