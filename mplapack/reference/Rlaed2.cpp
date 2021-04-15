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

void Rlaed2(INTEGER &k, INTEGER const n, INTEGER const n1, REAL *d, REAL *q, INTEGER const ldq, INTEGER *indxq, REAL &rho, REAL *z, REAL *dlamda, REAL *w, REAL *q2, INTEGER *indx, INTEGER *indxc, INTEGER *indxp, INTEGER *coltyp, INTEGER &info) {
    INTEGER n2 = 0;
    INTEGER n1p1 = 0;
    const REAL zero = 0.0;
    const REAL mone = -1.0;
    const REAL one = 1.0;
    const REAL two = 2.0;
    REAL t = 0.0;
    INTEGER i = 0;
    INTEGER imax = 0;
    INTEGER jmax = 0;
    REAL eps = 0.0;
    const REAL eight = 8.0;
    REAL tol = 0.0;
    INTEGER iq2 = 0;
    INTEGER j = 0;
    INTEGER k2 = 0;
    INTEGER nj = 0;
    INTEGER pj = 0;
    REAL s = 0.0;
    REAL c = 0.0;
    REAL tau = 0.0;
    INTEGER ctot[4];
    INTEGER ct = 0;
    INTEGER psm[4];
    INTEGER js = 0;
    INTEGER iq1 = 0;
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
    if (n < 0) {
        info = -2;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -6;
    } else if (min(1, (n / 2)) > n1 || (n / 2) < n1) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Rlaed2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    n2 = n - n1;
    n1p1 = n1 + 1;
    //
    if (rho < zero) {
        Rscal(n2, mone, &z[n1p1 - 1], 1);
    }
    //
    //     Normalize z so that norm(z) = 1.  Since z is the concatenation of
    //     two normalized vectors, norm2(z) = sqrt(2).
    //
    t = one / sqrt(two);
    Rscal(n, t, z, 1);
    //
    //     RHO = ABS( norm(z)**2 * RHO )
    //
    rho = abs(two * rho);
    //
    //     Sort the eigenvalues into increasing order
    //
    for (i = n1p1; i <= n; i = i + 1) {
        indxq[i - 1] += n1;
    }
    //
    //     re-integrate the deflated parts from the last pass
    //
    for (i = 1; i <= n; i = i + 1) {
        dlamda[i - 1] = d[indxq[i - 1] - 1];
    }
    Rlamrg(n1, n2, dlamda, 1, 1, indxc);
    for (i = 1; i <= n; i = i + 1) {
        indx[i - 1] = indxq[indxc[i - 1] - 1];
    }
    //
    //     Calculate the allowable deflation tolerance
    //
    imax = iRamax(n, z, 1);
    jmax = iRamax(n, d, 1);
    eps = Rlamch("Epsilon");
    tol = eight * eps * max(abs(d[jmax - 1]), abs(z[imax - 1]));
    //
    //     If the rank-1 modifier is small enough, no more needs to be done
    //     except to reorganize Q so that its columns correspond with the
    //     elements in D.
    //
    if (rho * abs(z[imax - 1]) <= tol) {
        k = 0;
        iq2 = 1;
        for (j = 1; j <= n; j = j + 1) {
            i = indx[j - 1];
            Rcopy(n, &q[(i - 1) * ldq], 1, &q2[iq2 - 1], 1);
            dlamda[j - 1] = d[i - 1];
            iq2 += n;
        }
        Rlacpy("A", n, n, q2, n, q, ldq);
        Rcopy(n, dlamda, 1, d, 1);
        goto statement_190;
    }
    //
    //     If there are multiple eigenvalues then the problem deflates.  Here
    //     the number of equal eigenvalues are found.  As each equal
    //     eigenvalue is found, an elementary reflector is computed to rotate
    //     the corresponding eigensubspace so that the corresponding
    //     components of Z are zero in this new basis.
    //
    for (i = 1; i <= n1; i = i + 1) {
        coltyp[i - 1] = 1;
    }
    for (i = n1p1; i <= n; i = i + 1) {
        coltyp[i - 1] = 3;
    }
    //
    k = 0;
    k2 = n + 1;
    for (j = 1; j <= n; j = j + 1) {
        nj = indx[j - 1];
        if (rho * abs(z[nj - 1]) <= tol) {
            //
            //           Deflate due to small z component.
            //
            k2 = k2 - 1;
            coltyp[nj - 1] = 4;
            indxp[k2 - 1] = nj;
            if (j == n) {
                goto statement_100;
            }
        } else {
            pj = nj;
            goto statement_80;
        }
    }
statement_80:
    j++;
    nj = indx[j - 1];
    if (j > n) {
        goto statement_100;
    }
    if (rho * abs(z[nj - 1]) <= tol) {
        //
        //        Deflate due to small z component.
        //
        k2 = k2 - 1;
        coltyp[nj - 1] = 4;
        indxp[k2 - 1] = nj;
    } else {
        //
        //        Check if eigenvalues are close enough to allow deflation.
        //
        s = z[pj - 1];
        c = z[nj - 1];
        //
        //        Find sqrt(a**2+b**2) without overflow or
        //        destructive underflow.
        //
        tau = Rlapy2(c, s);
        t = d[nj - 1] - d[pj - 1];
        c = c / tau;
        s = -s / tau;
        if (abs(t * c * s) <= tol) {
            //
            //           Deflation is possible.
            //
            z[nj - 1] = tau;
            z[pj - 1] = zero;
            if (coltyp[nj - 1] != coltyp[pj - 1]) {
                coltyp[nj - 1] = 2;
            }
            coltyp[pj - 1] = 4;
            Rrot(n, &q[(pj - 1) * ldq], 1, &q[(nj - 1) * ldq], 1, c, s);
            t = d[pj - 1] * pow2(c) + d[nj - 1] * pow2(s);
            d[nj - 1] = d[pj - 1] * pow2(s) + d[nj - 1] * pow2(c);
            d[pj - 1] = t;
            k2 = k2 - 1;
            i = 1;
        statement_90:
            if (k2 + i <= n) {
                if (d[pj - 1] < d[(indxp[(k2 + i) - 1]) - 1]) {
                    indxp[(k2 + i - 1) - 1] = indxp[(k2 + i) - 1];
                    indxp[(k2 + i) - 1] = pj;
                    i++;
                    goto statement_90;
                } else {
                    indxp[(k2 + i - 1) - 1] = pj;
                }
            } else {
                indxp[(k2 + i - 1) - 1] = pj;
            }
            pj = nj;
        } else {
            k++;
            dlamda[k - 1] = d[pj - 1];
            w[k - 1] = z[pj - 1];
            indxp[k - 1] = pj;
            pj = nj;
        }
    }
    goto statement_80;
statement_100:
    //
    //     Record the last eigenvalue.
    //
    k++;
    dlamda[k - 1] = d[pj - 1];
    w[k - 1] = z[pj - 1];
    indxp[k - 1] = pj;
    //
    //     Count up the total number of the various types of columns, then
    //     form a permutation which positions the four column types into
    //     four uniform groups (although one or more of these groups may be
    //     empty).
    //
    for (j = 1; j <= 4; j = j + 1) {
        ctot[j - 1] = 0;
    }
    for (j = 1; j <= n; j = j + 1) {
        ct = coltyp[j - 1];
        ctot[ct - 1]++;
    }
    //
    //     PSM(*) = Position in SubMatrix (of types 1 through 4)
    //
    psm[1 - 1] = 1;
    psm[2 - 1] = 1 + ctot[1 - 1];
    psm[3 - 1] = psm[2 - 1] + ctot[2 - 1];
    psm[4 - 1] = psm[3 - 1] + ctot[3 - 1];
    k = n - ctot[4 - 1];
    //
    //     Fill out the INDXC array so that the permutation which it induces
    //     will place all type-1 columns first, all type-2 columns next,
    //     then all type-3's, and finally all type-4's.
    //
    for (j = 1; j <= n; j = j + 1) {
        js = indxp[j - 1];
        ct = coltyp[js - 1];
        indx[psm[ct - 1] - 1] = js;
        indxc[psm[ct - 1] - 1] = j;
        psm[ct - 1]++;
    }
    //
    //     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
    //     and Q2 respectively.  The eigenvalues/vectors which were not
    //     deflated go into the first K slots of DLAMDA and Q2 respectively,
    //     while those which were deflated go into the last N - K slots.
    //
    i = 1;
    iq1 = 1;
    iq2 = 1 + (ctot[1 - 1] + ctot[2 - 1]) * n1;
    for (j = 1; j <= ctot[1 - 1]; j = j + 1) {
        js = indx[i - 1];
        Rcopy(n1, &q[(js - 1) * ldq], 1, &q2[iq1 - 1], 1);
        z[i - 1] = d[js - 1];
        i++;
        iq1 += n1;
    }
    //
    for (j = 1; j <= ctot[2 - 1]; j = j + 1) {
        js = indx[i - 1];
        Rcopy(n1, &q[(js - 1) * ldq], 1, &q2[iq1 - 1], 1);
        Rcopy(n2, &q[((n1 + 1) - 1) + (js - 1) * ldq], 1, &q2[iq2 - 1], 1);
        z[i - 1] = d[js - 1];
        i++;
        iq1 += n1;
        iq2 += n2;
    }
    //
    for (j = 1; j <= ctot[3 - 1]; j = j + 1) {
        js = indx[i - 1];
        Rcopy(n2, &q[((n1 + 1) - 1) + (js - 1) * ldq], 1, &q2[iq2 - 1], 1);
        z[i - 1] = d[js - 1];
        i++;
        iq2 += n2;
    }
    //
    iq1 = iq2;
    for (j = 1; j <= ctot[4 - 1]; j = j + 1) {
        js = indx[i - 1];
        Rcopy(n, &q[(js - 1) * ldq], 1, &q2[iq2 - 1], 1);
        iq2 += n;
        z[i - 1] = d[js - 1];
        i++;
    }
    //
    //     The deflated eigenvalues and their corresponding vectors go back
    //     into the last N - K slots of D and Q respectively.
    //
    if (k < n) {
        Rlacpy("A", n, ctot[4 - 1], &q2[iq1 - 1], n, &q[((k + 1) - 1) * ldq], ldq);
        Rcopy(n - k, &z[(k + 1) - 1], 1, &d[(k + 1) - 1], 1);
    }
    //
    //     Copy CTOT into COLTYP for referencing in Rlaed3.
    //
    for (j = 1; j <= 4; j = j + 1) {
        coltyp[j - 1] = ctot[j - 1];
    }
//
statement_190:;
    //
    //     End of Rlaed2
    //
}
