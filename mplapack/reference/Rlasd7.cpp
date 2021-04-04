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

void Rlasd7(INTEGER const &icompq, INTEGER const &nl, INTEGER const &nr, INTEGER const &sqre, INTEGER &k, REAL *d, REAL *z, REAL *zw, REAL *vf, REAL *vfw, REAL *vl, REAL *vlw, REAL const &alpha, REAL const &beta, REAL *dsigma, INTEGER *idx, arr_ref<INTEGER> idxp, arr_ref<INTEGER> idxq, arr_ref<INTEGER> perm, INTEGER &givptr, arr_ref<INTEGER, 2> givcol, INTEGER const &ldgcol, REAL *givnum, INTEGER const &ldgnum, REAL &c, REAL &s, INTEGER &info) {
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER nlp1 = 0;
    INTEGER nlp2 = 0;
    REAL z1 = 0.0;
    const REAL zero = 0.0;
    REAL tau = 0.0;
    INTEGER i = 0;
    INTEGER idxi = 0;
    REAL eps = 0.0;
    REAL tol = 0.0;
    const REAL eight = 8.0e+0;
    INTEGER k2 = 0;
    INTEGER j = 0;
    INTEGER jprev = 0;
    INTEGER idxjp = 0;
    INTEGER idxj = 0;
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
    //     .. Local Scalars ..
    //
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    n = nl + nr + 1;
    m = n + sqre;
    //
    if ((icompq < 0) || (icompq > 1)) {
        info = -1;
    } else if (nl < 1) {
        info = -2;
    } else if (nr < 1) {
        info = -3;
    } else if ((sqre < 0) || (sqre > 1)) {
        info = -4;
    } else if (ldgcol < n) {
        info = -22;
    } else if (ldgnum < n) {
        info = -24;
    }
    if (info != 0) {
        Mxerbla("Rlasd7", -info);
        return;
    }
    //
    nlp1 = nl + 1;
    nlp2 = nl + 2;
    if (icompq == 1) {
        givptr = 0;
    }
    //
    //     Generate the first part of the vector Z and move the singular
    //     values in the first part of D one position backward.
    //
    z1 = alpha * vl[nlp1 - 1];
    vl[nlp1 - 1] = zero;
    tau = vf[nlp1 - 1];
    for (i = nl; i >= 1; i = i - 1) {
        z[(i + 1) - 1] = alpha * vl[i - 1];
        vl[i - 1] = zero;
        vf[(i + 1) - 1] = vf[i - 1];
        d[(i + 1) - 1] = d[i - 1];
        idxq[(i + 1) - 1] = idxq[i - 1] + 1;
    }
    vf[1 - 1] = tau;
    //
    //     Generate the second part of the vector Z.
    //
    for (i = nlp2; i <= m; i = i + 1) {
        z[i - 1] = beta * vf[i - 1];
        vf[i - 1] = zero;
    }
    //
    //     Sort the singular values INTEGERo increasing order
    //
    for (i = nlp2; i <= n; i = i + 1) {
        idxq[i - 1] += nlp1;
    }
    //
    //     DSIGMA, IDXC, IDXC, and ZW are used as storage space.
    //
    for (i = 2; i <= n; i = i + 1) {
        dsigma[i - 1] = d[idxq[i - 1] - 1];
        zw[i - 1] = z[idxq[i - 1] - 1];
        vfw[i - 1] = vf[idxq[i - 1] - 1];
        vlw[i - 1] = vl[idxq[i - 1] - 1];
    }
    //
    Rlamrg(nl, nr, dsigma[2 - 1], 1, 1, idx[2 - 1]);
    //
    for (i = 2; i <= n; i = i + 1) {
        idxi = 1 + idx[i - 1];
        d[i - 1] = dsigma[idxi - 1];
        z[i - 1] = zw[idxi - 1];
        vf[i - 1] = vfw[idxi - 1];
        vl[i - 1] = vlw[idxi - 1];
    }
    //
    //     Calculate the allowable deflation tolerance
    //
    eps = dlamch("Epsilon");
    tol = max(abs(alpha), abs(beta));
    tol = eight * eight * eps * max(abs(d[n - 1]), tol);
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
            if (j == n) {
                goto statement_100;
            }
        } else {
            jprev = j;
            goto statement_70;
        }
    }
statement_70:
    j = jprev;
statement_80:
    j++;
    if (j > n) {
        goto statement_90;
    }
    if (abs(z[j - 1]) <= tol) {
        //
        //        Deflate due to small z component.
        //
        k2 = k2 - 1;
        idxp[k2 - 1] = j;
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
            tau = Rlapy2[(c - 1) + (s - 1) * ldRlapy2];
            z[j - 1] = tau;
            z[jprev - 1] = zero;
            c = c / tau;
            s = -s / tau;
            //
            //           Record the appropriate Givens rotation
            //
            if (icompq == 1) {
                givptr++;
                idxjp = idxq[(idx[jprev - 1] + 1) - 1];
                idxj = idxq[(idx[j - 1] + 1) - 1];
                if (idxjp <= nlp1) {
                    idxjp = idxjp - 1;
                }
                if (idxj <= nlp1) {
                    idxj = idxj - 1;
                }
                givcol[(givptr - 1) + (2 - 1) * ldgivcol] = idxjp;
                givcol[(givptr - 1)] = idxj;
                givnum[(givptr - 1) + (2 - 1) * ldgivnum] = c;
                givnum[(givptr - 1)] = s;
            }
            Rrot(1, vf[jprev - 1], 1, vf[j - 1], 1, c, s);
            Rrot(1, vl[jprev - 1], 1, vl[j - 1], 1, c, s);
            k2 = k2 - 1;
            idxp[k2 - 1] = jprev;
            jprev = j;
        } else {
            k++;
            zw[k - 1] = z[jprev - 1];
            dsigma[k - 1] = d[jprev - 1];
            idxp[k - 1] = jprev;
            jprev = j;
        }
    }
    goto statement_80;
statement_90:
    //
    //     Record the last singular value.
    //
    k++;
    zw[k - 1] = z[jprev - 1];
    dsigma[k - 1] = d[jprev - 1];
    idxp[k - 1] = jprev;
//
statement_100:
    //
    //     Sort the singular values INTEGERo DSIGMA. The singular values which
    //     were not deflated go INTEGERo the first K slots of DSIGMA, except
    //     that DSIGMA(1) is treated separately.
    //
    for (j = 2; j <= n; j = j + 1) {
        jp = idxp[j - 1];
        dsigma[j - 1] = d[jp - 1];
        vfw[j - 1] = vf[jp - 1];
        vlw[j - 1] = vl[jp - 1];
    }
    if (icompq == 1) {
        for (j = 2; j <= n; j = j + 1) {
            jp = idxp[j - 1];
            perm[j - 1] = idxq[(idx[jp - 1] + 1) - 1];
            if (perm[j - 1] <= nlp1) {
                perm[j - 1] = perm[j - 1] - 1;
            }
        }
    }
    //
    //     The deflated singular values go back INTEGERo the last N - K slots of
    //     D.
    //
    Rcopy(n - k, dsigma[(k + 1) - 1], 1, d[(k + 1) - 1], 1);
    //
    //     Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
    //     VL(M).
    //
    dsigma[1 - 1] = zero;
    hlftol = tol / two;
    if (abs(dsigma[2 - 1]) <= hlftol) {
        dsigma[2 - 1] = hlftol;
    }
    if (m > n) {
        z[1 - 1] = Rlapy2[(z1 - 1) + (z[m - 1] - 1) * ldRlapy2];
        if (z[1 - 1] <= tol) {
            c = one;
            s = zero;
            z[1 - 1] = tol;
        } else {
            c = z1 / z[1 - 1];
            s = -z[m - 1] / z[1 - 1];
        }
        Rrot(1, vf[m - 1], 1, vf[1 - 1], 1, c, s);
        Rrot(1, vl[m - 1], 1, vl[1 - 1], 1, c, s);
    } else {
        if (abs(z1) <= tol) {
            z[1 - 1] = tol;
        } else {
            z[1 - 1] = z1;
        }
    }
    //
    //     Restore Z, VF, and VL.
    //
    Rcopy(k - 1, zw[2 - 1], 1, z[2 - 1], 1);
    Rcopy(n - 1, vfw[2 - 1], 1, vf[2 - 1], 1);
    Rcopy(n - 1, vlw[2 - 1], 1, vl[2 - 1], 1);
    //
    //     End of Rlasd7
    //
}
