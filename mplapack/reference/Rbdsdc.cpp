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

void Rbdsdc(const char *uplo, const char *compq, INTEGER const &n, REAL *d, REAL *e, REAL *u, INTEGER const &ldu, REAL *vt, INTEGER const &ldvt, REAL *q, arr_ref<INTEGER> iq, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER iuplo = 0;
    INTEGER icompq = 0;
    INTEGER smlsiz = 0;
    const REAL one = 1.0;
    INTEGER nm1 = 0;
    INTEGER wstart = 0;
    INTEGER qstart = 0;
    INTEGER i = 0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL r = 0.0;
    const REAL zero = 0.0;
    INTEGER iu = 0;
    INTEGER ivt = 0;
    REAL orgnrm = 0.0;
    INTEGER ierr = 0;
    REAL eps = 0.0;
    const REAL two = 2.0e+0;
    INTEGER mlvl = 0;
    INTEGER smlszp = 0;
    INTEGER difl = 0;
    INTEGER difr = 0;
    INTEGER z = 0;
    INTEGER ic = 0;
    INTEGER is = 0;
    INTEGER poles = 0;
    INTEGER givnum = 0;
    INTEGER k = 0;
    INTEGER givptr = 0;
    INTEGER perm = 0;
    INTEGER givcol = 0;
    INTEGER start = 0;
    INTEGER sqre = 0;
    INTEGER nsize = 0;
    INTEGER ii = 0;
    INTEGER kk = 0;
    REAL p = 0.0;
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
    //  (N-1).  Sven, 17 Feb 05.
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
    iuplo = 0;
    if (Mlsame(uplo, "U")) {
        iuplo = 1;
    }
    if (Mlsame(uplo, "L")) {
        iuplo = 2;
    }
    if (Mlsame(compq, "N")) {
        icompq = 0;
    } else if (Mlsame(compq, "P")) {
        icompq = 1;
    } else if (Mlsame(compq, "I")) {
        icompq = 2;
    } else {
        icompq = -1;
    }
    if (iuplo == 0) {
        info = -1;
    } else if (icompq < 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if ((ldu < 1) || ((icompq == 2) && (ldu < n))) {
        info = -7;
    } else if ((ldvt < 1) || ((icompq == 2) && (ldvt < n))) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rbdsdc", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    smlsiz = iMlaenv[(9 - 1) + ("Rbdsdc" - 1) * ldiMlaenv];
    if (n == 1) {
        if (icompq == 1) {
            q[1 - 1] = sign[(one - 1) + (d[1 - 1] - 1) * ldsign];
            q[(1 + smlsiz * n) - 1] = one;
        } else if (icompq == 2) {
            u[(1 - 1)] = sign[(one - 1) + (d[1 - 1] - 1) * ldsign];
            vt[(1 - 1)] = one;
        }
        d[1 - 1] = abs(d[1 - 1]);
        return;
    }
    nm1 = n - 1;
    //
    //     If matrix lower bidiagonal, rotate to be upper bidiagonal
    //     by applying Givens rotations on the left
    //
    wstart = 1;
    qstart = 3;
    if (icompq == 1) {
        Rcopy(n, d, 1, q[1 - 1], 1);
        Rcopy(n - 1, e, 1, q[(n + 1) - 1], 1);
    }
    if (iuplo == 2) {
        qstart = 5;
        if (icompq == 2) {
            wstart = 2 * n - 1;
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            Rlartg(d[i - 1], e[i - 1], cs, sn, r);
            d[i - 1] = r;
            e[i - 1] = sn * d[(i + 1) - 1];
            d[(i + 1) - 1] = cs * d[(i + 1) - 1];
            if (icompq == 1) {
                q[(i + 2 * n) - 1] = cs;
                q[(i + 3 * n) - 1] = sn;
            } else if (icompq == 2) {
                work[i - 1] = cs;
                work[(nm1 + i) - 1] = -sn;
            }
        }
    }
    //
    //     If ICOMPQ = 0, use Rlasdq to compute the singular values.
    //
    if (icompq == 0) {
        //        Ignore WSTART, instead using WORK( 1 ), since the two vectors
        //        for CS and -SN above are added only if ICOMPQ == 2,
        //        and adding them exceeds documented WORK size of 4*n.
        Rlasdq("U", 0, n, 0, 0, 0, d, e, vt, ldvt, u, ldu, u, ldu, work[1 - 1], info);
        goto statement_40;
    }
    //
    //     If N is smaller than the minimum divide size SMLSIZ, then solve
    //     the problem with another solver.
    //
    if (n <= smlsiz) {
        if (icompq == 2) {
            Rlaset("A", n, n, zero, one, u, ldu);
            Rlaset("A", n, n, zero, one, vt, ldvt);
            Rlasdq("U", 0, n, n, n, 0, d, e, vt, ldvt, u, ldu, u, ldu, work[wstart - 1], info);
        } else if (icompq == 1) {
            iu = 1;
            ivt = iu + n;
            Rlaset("A", n, n, zero, one, q[(iu + (qstart - 1) * n) - 1], n);
            Rlaset("A", n, n, zero, one, q[(ivt + (qstart - 1) * n) - 1], n);
            Rlasdq("U", 0, n, n, n, 0, d, e, q[(ivt + (qstart - 1) * n) - 1], n, q[(iu + (qstart - 1) * n) - 1], n, q[(iu + (qstart - 1) * n) - 1], n, work[wstart - 1], info);
        }
        goto statement_40;
    }
    //
    if (icompq == 2) {
        Rlaset("A", n, n, zero, one, u, ldu);
        Rlaset("A", n, n, zero, one, vt, ldvt);
    }
    //
    //     Scale.
    //
    orgnrm = Rlanst[("M" - 1) + (n - 1) * ldRlanst];
    if (orgnrm == zero) {
        return;
    }
    Rlascl("G", 0, 0, orgnrm, one, n, 1, d, n, ierr);
    Rlascl("G", 0, 0, orgnrm, one, nm1, 1, e, nm1, ierr);
    //
    eps = (0.9e+0) * dlamch("Epsilon");
    //
    mlvl = INTEGER(log[(n.real() / smlsiz + 1.real()) - 1] / log[two - 1]) + 1;
    smlszp = smlsiz + 1;
    //
    if (icompq == 1) {
        iu = 1;
        ivt = 1 + smlsiz;
        difl = ivt + smlszp;
        difr = difl + mlvl;
        z = difr + mlvl * 2;
        ic = z + mlvl;
        is = ic + 1;
        poles = is + 1;
        givnum = poles + 2 * mlvl;
        //
        k = 1;
        givptr = 2;
        perm = 3;
        givcol = perm + mlvl;
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        if (abs(d[i - 1]) < eps) {
            d[i - 1] = sign[(eps - 1) + (d[i - 1] - 1) * ldsign];
        }
    }
    //
    start = 1;
    sqre = 0;
    //
    for (i = 1; i <= nm1; i = i + 1) {
        if ((abs(e[i - 1]) < eps) || (i == nm1)) {
            //
            //           Subproblem found. First determine its size and then
            //           apply divide and conquer on it.
            //
            if (i < nm1) {
                //
                //              A subproblem with E(I) small for I < NM1.
                //
                nsize = i - start + 1;
            } else if (abs(e[i - 1]) >= eps) {
                //
                //              A subproblem with E(NM1) not too small but I = NM1.
                //
                nsize = n - start + 1;
            } else {
                //
                //              A subproblem with E(NM1) small. This implies an
                //              1-by-1 subproblem at D(N). Solve this 1-by-1 problem
                //              first.
                //
                nsize = i - start + 1;
                if (icompq == 2) {
                    u[(n - 1) + (n - 1) * ldu] = sign[(one - 1) + (d[n - 1] - 1) * ldsign];
                    vt[(n - 1) + (n - 1) * ldvt] = one;
                } else if (icompq == 1) {
                    q[(n + (qstart - 1) * n) - 1] = sign[(one - 1) + (d[n - 1] - 1) * ldsign];
                    q[(n + (smlsiz + qstart - 1) * n) - 1] = one;
                }
                d[n - 1] = abs(d[n - 1]);
            }
            if (icompq == 2) {
                Rlasd0(nsize, sqre, d[start - 1], e[start - 1], u[(start - 1) + (start - 1) * ldu], ldu, vt[(start - 1) + (start - 1) * ldvt], ldvt, smlsiz, iwork, work[wstart - 1], info);
            } else {
                Rlasda(icompq, smlsiz, nsize, sqre, d[start - 1], e[start - 1], q[(start + (iu + qstart - 2) * n) - 1], n, q[(start + (ivt + qstart - 2) * n) - 1], iq[(start + k * n) - 1], q[(start + (difl + qstart - 2) * n) - 1], q[(start + (difr + qstart - 2) * n) - 1], q[(start + (z + qstart - 2) * n) - 1], q[(start + (poles + qstart - 2) * n) - 1], iq[(start + givptr * n) - 1], iq[(start + givcol * n) - 1], n, iq[(start + perm * n) - 1], q[(start + (givnum + qstart - 2) * n) - 1], q[(start + (ic + qstart - 2) * n) - 1], q[(start + (is + qstart - 2) * n) - 1], work[wstart - 1], iwork, info);
            }
            if (info != 0) {
                return;
            }
            start = i + 1;
        }
    }
    //
    //     Unscale
    //
    Rlascl("G", 0, 0, one, orgnrm, n, 1, d, n, ierr);
statement_40:
    //
    //     Use Selection Sort to minimize swaps of singular vectors
    //
    for (ii = 2; ii <= n; ii = ii + 1) {
        i = ii - 1;
        kk = i;
        p = d[i - 1];
        for (j = ii; j <= n; j = j + 1) {
            if (d[j - 1] > p) {
                kk = j;
                p = d[j - 1];
            }
        }
        if (kk != i) {
            d[kk - 1] = d[i - 1];
            d[i - 1] = p;
            if (icompq == 1) {
                iq[i - 1] = kk;
            } else if (icompq == 2) {
                Rswap(n, u[(i - 1) * ldu], 1, u[(kk - 1) * ldu], 1);
                Rswap(n, vt[(i - 1)], ldvt, vt[(kk - 1)], ldvt);
            }
        } else if (icompq == 1) {
            iq[i - 1] = i;
        }
    }
    //
    //     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
    //
    if (icompq == 1) {
        if (iuplo == 1) {
            iq[n - 1] = 1;
        } else {
            iq[n - 1] = 0;
        }
    }
    //
    //     If B is lower bidiagonal, update U by those Givens rotations
    //     which rotated B to be upper bidiagonal
    //
    if ((iuplo == 2) && (icompq == 2)) {
        Rlasr("L", "V", "B", n, n, work[1 - 1], work[n - 1], u, ldu);
    }
    //
    //     End of Rbdsdc
    //
}
