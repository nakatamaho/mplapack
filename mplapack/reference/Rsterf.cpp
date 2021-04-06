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

void Rsterf(INTEGER const n, REAL *d, REAL *e, INTEGER &info) {
    REAL eps = 0.0;
    REAL eps2 = 0.0;
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    const REAL three = 3.0;
    REAL ssfmax = 0.0;
    REAL ssfmin = 0.0;
    REAL rmax = 0.0;
    const INTEGER maxit = 30;
    INTEGER nmaxit = 0;
    const REAL zero = 0.0;
    REAL sigma = 0.0;
    INTEGER jtot = 0;
    INTEGER l1 = 0;
    INTEGER m = 0;
    INTEGER l = 0;
    INTEGER lsv = 0;
    INTEGER lend = 0;
    INTEGER lendsv = 0;
    REAL anorm = 0.0;
    INTEGER iscale = 0;
    INTEGER i = 0;
    REAL p = 0.0;
    REAL rte = 0.0;
    REAL rt1 = 0.0;
    REAL rt2 = 0.0;
    const REAL two = 2.0;
    REAL r = 0.0;
    REAL c = 0.0;
    REAL s = 0.0;
    REAL gamma = 0.0;
    REAL bb = 0.0;
    REAL oldc = 0.0;
    REAL oldgam = 0.0;
    REAL alpha = 0.0;
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
    //     Quick return if possible
    //
    if (n < 0) {
        info = -1;
        Mxerbla("Rsterf", -info);
        return;
    }
    if (n <= 1) {
        return;
    }
    //
    //     Determine the unit roundoff for this environment.
    //
    eps = Rlamch("E");
    eps2 = pow2(eps);
    safmin = Rlamch("S");
    safmax = one / safmin;
    ssfmax = sqrt(safmax) / three;
    ssfmin = sqrt(safmin) / eps2;
    rmax = Rlamch("O");
    //
    //     Compute the eigenvalues of the tridiagonal matrix.
    //
    nmaxit = n * maxit;
    sigma = zero;
    jtot = 0;
    //
    //     Determine where the matrix splits and choose QL or QR iteration
    //     for each block, according to whether top or bottom diagonal
    //     element is smaller.
    //
    l1 = 1;
//
statement_10:
    if (l1 > n) {
        goto statement_170;
    }
    if (l1 > 1) {
        e[(l1 - 1) - 1] = zero;
    }
    for (m = l1; m <= n - 1; m = m + 1) {
        if (abs(e[m - 1]) <= (sqrt(abs(d[m - 1])) * sqrt(abs(d[(m + 1) - 1]))) * eps) {
            e[m - 1] = zero;
            goto statement_30;
        }
    }
    m = n;
//
statement_30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
        goto statement_10;
    }
    //
    //     Scale submatrix in rows and columns L to LEND
    //
    anorm = Rlanst("M", lend - l + 1, &d[l - 1], &e[l - 1]);
    iscale = 0;
    if (anorm == zero) {
        goto statement_10;
    }
    if ((anorm > ssfmax)) {
        iscale = 1;
        Rlascl("G", 0, 0, anorm, ssfmax, lend - l + 1, 1, &d[l - 1], n, info);
        Rlascl("G", 0, 0, anorm, ssfmax, lend - l, 1, &e[l - 1], n, info);
    } else if (anorm < ssfmin) {
        iscale = 2;
        Rlascl("G", 0, 0, anorm, ssfmin, lend - l + 1, 1, &d[l - 1], n, info);
        Rlascl("G", 0, 0, anorm, ssfmin, lend - l, 1, &e[l - 1], n, info);
    }
    //
    for (i = l; i <= lend - 1; i = i + 1) {
        e[i - 1] = pow2(e[i - 1]);
    }
    //
    //     Choose between QL and QR iteration
    //
    if (abs(d[lend - 1]) < abs(d[l - 1])) {
        lend = lsv;
        l = lendsv;
    }
    //
    if (lend >= l) {
    //
    //        QL Iteration
    //
    //        Look for small subdiagonal element.
    //
    statement_50:
        if (l != lend) {
            for (m = l; m <= lend - 1; m = m + 1) {
                if (abs(e[m - 1]) <= eps2 * abs(d[m - 1] * d[(m + 1) - 1])) {
                    goto statement_70;
                }
            }
        }
        m = lend;
    //
    statement_70:
        if (m < lend) {
            e[m - 1] = zero;
        }
        p = d[l - 1];
        if (m == l) {
            goto statement_90;
        }
        //
        //        If remaining matrix is 2 by 2, use Rlae2 to compute its
        //        eigenvalues.
        //
        if (m == l + 1) {
            rte = sqrt(e[l - 1]);
            Rlae2(d[l - 1], rte, d[(l + 1) - 1], rt1, rt2);
            d[l - 1] = rt1;
            d[(l + 1) - 1] = rt2;
            e[l - 1] = zero;
            l += 2;
            if (l <= lend) {
                goto statement_50;
            }
            goto statement_150;
        }
        //
        if (jtot == nmaxit) {
            goto statement_150;
        }
        jtot++;
        //
        //        Form shift.
        //
        rte = sqrt(e[l - 1]);
        sigma = (d[(l + 1) - 1] - p) / (two * rte);
        r = Rlapy2(sigma, one);
        sigma = p - (rte / (sigma + sign(r, sigma)));
        //
        c = one;
        s = zero;
        gamma = d[m - 1] - sigma;
        p = gamma * gamma;
        //
        //        Inner loop
        //
        for (i = m - 1; i >= l; i = i - 1) {
            bb = e[i - 1];
            r = p + bb;
            if (i != m - 1) {
                e[(i + 1) - 1] = s * r;
            }
            oldc = c;
            c = p / r;
            s = bb / r;
            oldgam = gamma;
            alpha = d[i - 1];
            gamma = c * (alpha - sigma) - s * oldgam;
            d[(i + 1) - 1] = oldgam + (alpha - gamma);
            if (c != zero) {
                p = (gamma * gamma) / c;
            } else {
                p = oldc * bb;
            }
        }
        //
        e[l - 1] = s * p;
        d[l - 1] = sigma + gamma;
        goto statement_50;
    //
    //        Eigenvalue found.
    //
    statement_90:
        d[l - 1] = p;
        //
        l++;
        if (l <= lend) {
            goto statement_50;
        }
        goto statement_150;
        //
    } else {
    //
    //        QR Iteration
    //
    //        Look for small superdiagonal element.
    //
    statement_100:
        for (m = l; m >= lend + 1; m = m - 1) {
            if (abs(e[(m - 1) - 1]) <= eps2 * abs(d[m - 1] * d[(m - 1) - 1])) {
                goto statement_120;
            }
        }
        m = lend;
    //
    statement_120:
        if (m > lend) {
            e[(m - 1) - 1] = zero;
        }
        p = d[l - 1];
        if (m == l) {
            goto statement_140;
        }
        //
        //        If remaining matrix is 2 by 2, use Rlae2 to compute its
        //        eigenvalues.
        //
        if (m == l - 1) {
            rte = sqrt(e[(l - 1) - 1]);
            Rlae2(d[l - 1], rte, d[(l - 1) - 1], rt1, rt2);
            d[l - 1] = rt1;
            d[(l - 1) - 1] = rt2;
            e[(l - 1) - 1] = zero;
            l = l - 2;
            if (l >= lend) {
                goto statement_100;
            }
            goto statement_150;
        }
        //
        if (jtot == nmaxit) {
            goto statement_150;
        }
        jtot++;
        //
        //        Form shift.
        //
        rte = sqrt(e[(l - 1) - 1]);
        sigma = (d[(l - 1) - 1] - p) / (two * rte);
        r = Rlapy2(sigma, one);
        sigma = p - (rte / (sigma + sign(r, sigma)));
        //
        c = one;
        s = zero;
        gamma = d[m - 1] - sigma;
        p = gamma * gamma;
        //
        //        Inner loop
        //
        for (i = m; i <= l - 1; i = i + 1) {
            bb = e[i - 1];
            r = p + bb;
            if (i != m) {
                e[(i - 1) - 1] = s * r;
            }
            oldc = c;
            c = p / r;
            s = bb / r;
            oldgam = gamma;
            alpha = d[(i + 1) - 1];
            gamma = c * (alpha - sigma) - s * oldgam;
            d[i - 1] = oldgam + (alpha - gamma);
            if (c != zero) {
                p = (gamma * gamma) / c;
            } else {
                p = oldc * bb;
            }
        }
        //
        e[(l - 1) - 1] = s * p;
        d[l - 1] = sigma + gamma;
        goto statement_100;
    //
    //        Eigenvalue found.
    //
    statement_140:
        d[l - 1] = p;
        //
        l = l - 1;
        if (l >= lend) {
            goto statement_100;
        }
        goto statement_150;
        //
    }
//
//     Undo scaling if necessary
//
statement_150:
    if (iscale == 1) {
        Rlascl("G", 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, &d[lsv - 1], n, info);
    }
    if (iscale == 2) {
        Rlascl("G", 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, &d[lsv - 1], n, info);
    }
    //
    //     Check for no convergence to an eigenvalue after a total
    //     of N*MAXIT iterations.
    //
    if (jtot < nmaxit) {
        goto statement_10;
    }
    for (i = 1; i <= n - 1; i = i + 1) {
        if (e[i - 1] != zero) {
            info++;
        }
    }
    goto statement_180;
//
//     Sort eigenvalues in increasing order.
//
statement_170:
    Rlasrt("I", n, d, info);
//
statement_180:;
    //
    //     End of Rsterf
    //
}
