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

void Rsteqr(const char *compz, INTEGER const &n, REAL *d, REAL *e, REAL *z, INTEGER const &ldz, REAL *work, INTEGER &info) {
    INTEGER icompz = 0;
    const REAL one = 1.0;
    REAL eps = 0.0;
    REAL eps2 = 0.0;
    REAL safmin = 0.0;
    REAL safmax = 0.0;
    const REAL three = 3.0;
    REAL ssfmax = 0.0;
    REAL ssfmin = 0.0;
    const REAL zero = 0.0;
    const INTEGER maxit = 30;
    INTEGER nmaxit = 0;
    INTEGER jtot = 0;
    INTEGER l1 = 0;
    INTEGER nm1 = 0;
    INTEGER m = 0;
    REAL tst = 0.0;
    INTEGER l = 0;
    INTEGER lsv = 0;
    INTEGER lend = 0;
    INTEGER lendsv = 0;
    REAL anorm = 0.0;
    INTEGER iscale = 0;
    INTEGER lendm1 = 0;
    REAL p = 0.0;
    REAL rt1 = 0.0;
    REAL rt2 = 0.0;
    REAL c = 0.0;
    REAL s = 0.0;
    const REAL two = 2.0;
    REAL g = 0.0;
    REAL r = 0.0;
    INTEGER mm1 = 0;
    INTEGER i = 0;
    REAL f = 0.0;
    REAL b = 0.0;
    INTEGER mm = 0;
    INTEGER lendp1 = 0;
    INTEGER lm1 = 0;
    INTEGER ii = 0;
    INTEGER k = 0;
    INTEGER j = 0;
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
    if (Mlsame(compz, "N")) {
        icompz = 0;
    } else if (Mlsame(compz, "V")) {
        icompz = 1;
    } else if (Mlsame(compz, "I")) {
        icompz = 2;
    } else {
        icompz = -1;
    }
    if (icompz < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if ((ldz < 1) || (icompz > 0 && ldz < max((INTEGER)1, n))) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rsteqr", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    if (n == 1) {
        if (icompz == 2) {
            z[(1 - 1)] = one;
        }
        return;
    }
    //
    //     Determine the unit roundoff and over/underflow thresholds.
    //
    eps = dlamch("E");
    eps2 = pow2(eps);
    safmin = dlamch("S");
    safmax = one / safmin;
    ssfmax = sqrt(safmax) / three;
    ssfmin = sqrt(safmin) / eps2;
    //
    //     Compute the eigenvalues and eigenvectors of the tridiagonal
    //     matrix.
    //
    if (icompz == 2) {
        Rlaset("Full", n, n, zero, one, z, ldz);
    }
    //
    nmaxit = n * maxit;
    jtot = 0;
    //
    //     Determine where the matrix splits and choose QL or QR iteration
    //     for each block, according to whether top or bottom diagonal
    //     element is smaller.
    //
    l1 = 1;
    nm1 = n - 1;
//
statement_10:
    if (l1 > n) {
        goto statement_160;
    }
    if (l1 > 1) {
        e[(l1 - 1) - 1] = zero;
    }
    if (l1 <= nm1) {
        for (m = l1; m <= nm1; m = m + 1) {
            tst = abs(e[m - 1]);
            if (tst == zero) {
                goto statement_30;
            }
            if (tst <= (sqrt(abs(d[m - 1])) * sqrt(abs(d[(m + 1) - 1]))) * eps) {
                e[m - 1] = zero;
                goto statement_30;
            }
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
    anorm = Rlanst[("M" - 1) + ((lend - l + 1) - 1) * ldRlanst];
    iscale = 0;
    if (anorm == zero) {
        goto statement_10;
    }
    if (anorm > ssfmax) {
        iscale = 1;
        Rlascl("G", 0, 0, anorm, ssfmax, lend - l + 1, 1, d[l - 1], n, info);
        Rlascl("G", 0, 0, anorm, ssfmax, lend - l, 1, e[l - 1], n, info);
    } else if (anorm < ssfmin) {
        iscale = 2;
        Rlascl("G", 0, 0, anorm, ssfmin, lend - l + 1, 1, d[l - 1], n, info);
        Rlascl("G", 0, 0, anorm, ssfmin, lend - l, 1, e[l - 1], n, info);
    }
    //
    //     Choose between QL and QR iteration
    //
    if (abs(d[lend - 1]) < abs(d[l - 1])) {
        lend = lsv;
        l = lendsv;
    }
    //
    if (lend > l) {
    //
    //        QL Iteration
    //
    //        Look for small subdiagonal element.
    //
    statement_40:
        if (l != lend) {
            lendm1 = lend - 1;
            for (m = l; m <= lendm1; m = m + 1) {
                tst = pow2(abs(e[m - 1]));
                if (tst <= (eps2 * abs(d[m - 1])) * abs(d[(m + 1) - 1]) + safmin) {
                    goto statement_60;
                }
            }
        }
        //
        m = lend;
    //
    statement_60:
        if (m < lend) {
            e[m - 1] = zero;
        }
        p = d[l - 1];
        if (m == l) {
            goto statement_80;
        }
        //
        //        If remaining matrix is 2-by-2, use Rlae2 or SLAEV2
        //        to compute its eigensystem.
        //
        if (m == l + 1) {
            if (icompz > 0) {
                Rlaev2(d[l - 1], e[l - 1], d[(l + 1) - 1], rt1, rt2, c, s);
                work[l - 1] = c;
                work[(n - 1 + l) - 1] = s;
                Rlasr("R", "V", "B", n, 2, work[l - 1], work[(n - 1 + l) - 1], z[(l - 1) * ldz], ldz);
            } else {
                Rlae2(d[l - 1], e[l - 1], d[(l + 1) - 1], rt1, rt2);
            }
            d[l - 1] = rt1;
            d[(l + 1) - 1] = rt2;
            e[l - 1] = zero;
            l += 2;
            if (l <= lend) {
                goto statement_40;
            }
            goto statement_140;
        }
        //
        if (jtot == nmaxit) {
            goto statement_140;
        }
        jtot++;
        //
        //        Form shift.
        //
        g = (d[(l + 1) - 1] - p) / (two * e[l - 1]);
        r = Rlapy2[(g - 1) + (one - 1) * ldRlapy2];
        g = d[m - 1] - p + (e[l - 1] / (g + sign[(r - 1) + (g - 1) * ldsign]));
        //
        s = one;
        c = one;
        p = zero;
        //
        //        Inner loop
        //
        mm1 = m - 1;
        for (i = mm1; i >= l; i = i - 1) {
            f = s * e[i - 1];
            b = c * e[i - 1];
            Rlartg(g, f, c, s, r);
            if (i != m - 1) {
                e[(i + 1) - 1] = r;
            }
            g = d[(i + 1) - 1] - p;
            r = (d[i - 1] - g) * s + two * c * b;
            p = s * r;
            d[(i + 1) - 1] = g + p;
            g = c * r - b;
            //
            //           If eigenvectors are desired, then save rotations.
            //
            if (icompz > 0) {
                work[i - 1] = c;
                work[(n - 1 + i) - 1] = -s;
            }
            //
        }
        //
        //        If eigenvectors are desired, then apply saved rotations.
        //
        if (icompz > 0) {
            mm = m - l + 1;
            Rlasr("R", "V", "B", n, mm, work[l - 1], work[(n - 1 + l) - 1], z[(l - 1) * ldz], ldz);
        }
        //
        d[l - 1] = d[l - 1] - p;
        e[l - 1] = g;
        goto statement_40;
    //
    //        Eigenvalue found.
    //
    statement_80:
        d[l - 1] = p;
        //
        l++;
        if (l <= lend) {
            goto statement_40;
        }
        goto statement_140;
        //
    } else {
    //
    //        QR Iteration
    //
    //        Look for small superdiagonal element.
    //
    statement_90:
        if (l != lend) {
            lendp1 = lend + 1;
            for (m = l; m >= lendp1; m = m - 1) {
                tst = pow2(abs(e[(m - 1) - 1]));
                if (tst <= (eps2 * abs(d[m - 1])) * abs(d[(m - 1) - 1]) + safmin) {
                    goto statement_110;
                }
            }
        }
        //
        m = lend;
    //
    statement_110:
        if (m > lend) {
            e[(m - 1) - 1] = zero;
        }
        p = d[l - 1];
        if (m == l) {
            goto statement_130;
        }
        //
        //        If remaining matrix is 2-by-2, use Rlae2 or SLAEV2
        //        to compute its eigensystem.
        //
        if (m == l - 1) {
            if (icompz > 0) {
                Rlaev2(d[(l - 1) - 1], e[(l - 1) - 1], d[l - 1], rt1, rt2, c, s);
                work[m - 1] = c;
                work[(n - 1 + m) - 1] = s;
                Rlasr("R", "V", "F", n, 2, work[m - 1], work[(n - 1 + m) - 1], z[((l - 1) - 1) * ldz], ldz);
            } else {
                Rlae2(d[(l - 1) - 1], e[(l - 1) - 1], d[l - 1], rt1, rt2);
            }
            d[(l - 1) - 1] = rt1;
            d[l - 1] = rt2;
            e[(l - 1) - 1] = zero;
            l = l - 2;
            if (l >= lend) {
                goto statement_90;
            }
            goto statement_140;
        }
        //
        if (jtot == nmaxit) {
            goto statement_140;
        }
        jtot++;
        //
        //        Form shift.
        //
        g = (d[(l - 1) - 1] - p) / (two * e[(l - 1) - 1]);
        r = Rlapy2[(g - 1) + (one - 1) * ldRlapy2];
        g = d[m - 1] - p + (e[(l - 1) - 1] / (g + sign[(r - 1) + (g - 1) * ldsign]));
        //
        s = one;
        c = one;
        p = zero;
        //
        //        Inner loop
        //
        lm1 = l - 1;
        for (i = m; i <= lm1; i = i + 1) {
            f = s * e[i - 1];
            b = c * e[i - 1];
            Rlartg(g, f, c, s, r);
            if (i != m) {
                e[(i - 1) - 1] = r;
            }
            g = d[i - 1] - p;
            r = (d[(i + 1) - 1] - g) * s + two * c * b;
            p = s * r;
            d[i - 1] = g + p;
            g = c * r - b;
            //
            //           If eigenvectors are desired, then save rotations.
            //
            if (icompz > 0) {
                work[i - 1] = c;
                work[(n - 1 + i) - 1] = s;
            }
            //
        }
        //
        //        If eigenvectors are desired, then apply saved rotations.
        //
        if (icompz > 0) {
            mm = l - m + 1;
            Rlasr("R", "V", "F", n, mm, work[m - 1], work[(n - 1 + m) - 1], z[(m - 1) * ldz], ldz);
        }
        //
        d[l - 1] = d[l - 1] - p;
        e[lm1 - 1] = g;
        goto statement_90;
    //
    //        Eigenvalue found.
    //
    statement_130:
        d[l - 1] = p;
        //
        l = l - 1;
        if (l >= lend) {
            goto statement_90;
        }
        goto statement_140;
        //
    }
//
//     Undo scaling if necessary
//
statement_140:
    if (iscale == 1) {
        Rlascl("G", 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d[lsv - 1], n, info);
        Rlascl("G", 0, 0, ssfmax, anorm, lendsv - lsv, 1, e[lsv - 1], n, info);
    } else if (iscale == 2) {
        Rlascl("G", 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d[lsv - 1], n, info);
        Rlascl("G", 0, 0, ssfmin, anorm, lendsv - lsv, 1, e[lsv - 1], n, info);
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
    goto statement_190;
//
//     Order eigenvalues and eigenvectors.
//
statement_160:
    if (icompz == 0) {
        //
        //        Use Quick Sort
        //
        Rlasrt("I", n, d, info);
        //
    } else {
        //
        //        Use Selection Sort to minimize swaps of eigenvectors
        //
        for (ii = 2; ii <= n; ii = ii + 1) {
            i = ii - 1;
            k = i;
            p = d[i - 1];
            for (j = ii; j <= n; j = j + 1) {
                if (d[j - 1] < p) {
                    k = j;
                    p = d[j - 1];
                }
            }
            if (k != i) {
                d[k - 1] = d[i - 1];
                d[i - 1] = p;
                Rswap(n, z[(i - 1) * ldz], 1, z[(k - 1) * ldz], 1);
            }
        }
    }
//
statement_190:;
    //
    //     End of Rsteqr
    //
}
