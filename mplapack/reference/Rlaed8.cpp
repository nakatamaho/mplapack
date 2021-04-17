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

void Rlaed8(INTEGER const icompq, INTEGER &k, INTEGER const n, INTEGER const qsiz, REAL *d, REAL *q, INTEGER const ldq, INTEGER *indxq, REAL &rho, INTEGER const cutpnt, REAL *z, REAL *dlamda, REAL *q2, INTEGER const ldq2, REAL *w, INTEGER *perm, INTEGER &givptr, INTEGER *givcol, REAL *givnum, INTEGER *indxp, INTEGER *indx, INTEGER &info) {
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER n1p1 = 0;
    const REAL zero = 0.0;
    const REAL mone = -1.0;
    const REAL one = 1.0;
    const REAL two = 2.0;
    REAL t = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER imax = 0;
    INTEGER jmax = 0;
    REAL eps = 0.0;
    const REAL eight = 8.0;
    REAL tol = 0.0;
    INTEGER k2 = 0;
    INTEGER jlam = 0;
    REAL s = 0.0;
    REAL c = 0.0;
    REAL tau = 0.0;
    INTEGER jp = 0;
    INTEGER ldgivcol = 2;
    INTEGER ldgivnum = 2;
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
    //
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
    if (icompq < 0 || icompq > 1) {
        info = -1;
    } else if (n < 0) {
        info = -3;
    } else if (icompq == 1 && qsiz < n) {
        info = -4;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -7;
    } else if (cutpnt < min((INTEGER)1, n) || cutpnt > n) {
        info = -10;
    } else if (ldq2 < max((INTEGER)1, n)) {
        info = -14;
    }
    if (info != 0) {
        Mxerbla("Rlaed8", -info);
        return;
    }
    //
    //     Need to initialize GIVPTR to O here in case of quick exit
    //     to prevent an unspecified code behavior (usually sigfault)
    //     when IWORK array on entry to *stedc is not zeroed
    //     (or at least some IWORK entries which used in *laed7 for GIVPTR).
    //
    givptr = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    n1 = cutpnt;
    n2 = n - n1;
    n1p1 = n1 + 1;
    //
    if (rho < zero) {
        Rscal(n2, mone, &z[n1p1 - 1], 1);
    }
    //
    //     Normalize z so that norm(z) = 1
    //
    t = one / sqrt(two);
    for (j = 1; j <= n; j = j + 1) {
        indx[j - 1] = j;
    }
    Rscal(n, t, z, 1);
    rho = abs(two * rho);
    //
    //     Sort the eigenvalues into increasing order
    //
    for (i = cutpnt + 1; i <= n; i = i + 1) {
        indxq[i - 1] += cutpnt;
    }
    for (i = 1; i <= n; i = i + 1) {
        dlamda[i - 1] = d[indxq[i - 1] - 1];
        w[i - 1] = z[indxq[i - 1] - 1];
    }
    i = 1;
    j = cutpnt + 1;
    Rlamrg(n1, n2, dlamda, 1, 1, indx);
    for (i = 1; i <= n; i = i + 1) {
        d[i - 1] = dlamda[indx[i - 1] - 1];
        z[i - 1] = w[indx[i - 1] - 1];
    }
    //
    //     Calculate the allowable deflation tolerance
    //
    imax = iRamax(n, z, 1);
    jmax = iRamax(n, d, 1);
    eps = Rlamch("Epsilon");
    tol = eight * eps * abs(d[jmax - 1]);
    //
    //     If the rank-1 modifier is small enough, no more needs to be done
    //     except to reorganize Q so that its columns correspond with the
    //     elements in D.
    //
    if (rho * abs(z[imax - 1]) <= tol) {
        k = 0;
        if (icompq == 0) {
            for (j = 1; j <= n; j = j + 1) {
                perm[j - 1] = indxq[indx[j - 1] - 1];
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                perm[j - 1] = indxq[indx[j - 1] - 1];
                Rcopy(qsiz, &q[(perm[j - 1] - 1) * ldq], 1, &q2[(j - 1) * ldq2], 1);
            }
            Rlacpy("A", qsiz, n, &q2[(1 - 1)], ldq2, &q[(1 - 1)], ldq);
        }
        return;
    }
    //
    //     If there are multiple eigenvalues then the problem deflates.  Here
    //     the number of equal eigenvalues are found.  As each equal
    //     eigenvalue is found, an elementary reflector is computed to rotate
    //     the corresponding eigensubspace so that the corresponding
    //     components of Z are zero in this new basis.
    //
    k = 0;
    k2 = n + 1;
    for (j = 1; j <= n; j = j + 1) {
        if (rho * abs(z[j - 1]) <= tol) {
            //
            //           Deflate due to small z component.
            //
            k2 = k2 - 1;
            indxp[k2 - 1] = j;
            if (j == n) {
                goto statement_110;
            }
        } else {
            jlam = j;
            goto statement_80;
        }
    }
statement_80:
    j++;
    if (j > n) {
        goto statement_100;
    }
    if (rho * abs(z[j - 1]) <= tol) {
        //
        //        Deflate due to small z component.
        //
        k2 = k2 - 1;
        indxp[k2 - 1] = j;
    } else {
        //
        //        Check if eigenvalues are close enough to allow deflation.
        //
        s = z[jlam - 1];
        c = z[j - 1];
        //
        //        Find sqrt(a**2+b**2) without overflow or
        //        destructive underflow.
        //
        tau = Rlapy2(c, s);
        t = d[j - 1] - d[jlam - 1];
        c = c / tau;
        s = -s / tau;
        if (abs(t * c * s) <= tol) {
            //
            //           Deflation is possible.
            //
            z[j - 1] = tau;
            z[jlam - 1] = zero;
            //
            //           Record the appropriate Givens rotation
            //
            givptr++;
            givcol[(givptr - 1) * ldgivcol] = indxq[indx[jlam - 1] - 1];
            givcol[(2 - 1) + (givptr - 1) * ldgivcol] = indxq[indx[j - 1] - 1];
            givnum[(givptr - 1) * ldgivnum] = c;
            givnum[(2 - 1) + (givptr - 1) * ldgivnum] = s;
            if (icompq == 1) {
                Rrot(qsiz, &q[(indxq[indx[jlam - 1] - 1] - 1) * ldq], 1, &q[(indxq[indx[j - 1] - 1] - 1) * ldq], 1, c, s);
            }
            t = d[jlam - 1] * c * c + d[j - 1] * s * s;
            d[j - 1] = d[jlam - 1] * s * s + d[j - 1] * c * c;
            d[jlam - 1] = t;
            k2 = k2 - 1;
            i = 1;
        statement_90:
            if (k2 + i <= n) {
                if (d[jlam - 1] < d[(indxp[(k2 + i) - 1]) - 1]) {
                    indxp[(k2 + i - 1) - 1] = indxp[(k2 + i) - 1];
                    indxp[(k2 + i) - 1] = jlam;
                    i++;
                    goto statement_90;
                } else {
                    indxp[(k2 + i - 1) - 1] = jlam;
                }
            } else {
                indxp[(k2 + i - 1) - 1] = jlam;
            }
            jlam = j;
        } else {
            k++;
            w[k - 1] = z[jlam - 1];
            dlamda[k - 1] = d[jlam - 1];
            indxp[k - 1] = jlam;
            jlam = j;
        }
    }
    goto statement_80;
statement_100:
    //
    //     Record the last eigenvalue.
    //
    k++;
    w[k - 1] = z[jlam - 1];
    dlamda[k - 1] = d[jlam - 1];
    indxp[k - 1] = jlam;
//
statement_110:
    //
    //     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
    //     and Q2 respectively.  The eigenvalues/vectors which were not
    //     deflated go into the first K slots of DLAMDA and Q2 respectively,
    //     while those which were deflated go into the last N - K slots.
    //
    if (icompq == 0) {
        for (j = 1; j <= n; j = j + 1) {
            jp = indxp[j - 1];
            dlamda[j - 1] = d[jp - 1];
            perm[j - 1] = indxq[indx[jp - 1] - 1];
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            jp = indxp[j - 1];
            dlamda[j - 1] = d[jp - 1];
            perm[j - 1] = indxq[indx[jp - 1] - 1];
            Rcopy(qsiz, &q[(perm[j - 1] - 1) * ldq], 1, &q2[(j - 1) * ldq2], 1);
        }
    }
    //
    //     The deflated eigenvalues and their corresponding vectors go back
    //     into the last N - K slots of D and Q respectively.
    //
    if (k < n) {
        if (icompq == 0) {
            Rcopy(n - k, &dlamda[(k + 1) - 1], 1, &d[(k + 1) - 1], 1);
        } else {
            Rcopy(n - k, &dlamda[(k + 1) - 1], 1, &d[(k + 1) - 1], 1);
            Rlacpy("A", qsiz, n - k, &q2[((k + 1) - 1) * ldq2], ldq2, &q[((k + 1) - 1) * ldq], ldq);
        }
    }
    //
    //     End of Rlaed8
    //
}
