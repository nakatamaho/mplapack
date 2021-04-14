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

inline REAL Rlamc3(REAL a, REAL b) {
    REAL c = a + b;
    return c;
}

void Rlasd3(INTEGER const nl, INTEGER const nr, INTEGER const sqre, INTEGER const k, REAL *d, REAL *q, INTEGER const ldq, REAL *dsigma, REAL *u, INTEGER const ldu, REAL *u2, INTEGER const ldu2, REAL *vt, INTEGER const ldvt, REAL *vt2, INTEGER const ldvt2, INTEGER *idxc, INTEGER *ctot, REAL *z, INTEGER &info) {
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER nlp1 = 0;
    INTEGER nlp2 = 0;
    const REAL zero = 0.0;
    INTEGER i = 0;
    REAL rho = 0.0;
    const REAL one = 1.0;
    INTEGER j = 0;
    const REAL negone = -1.0;
    REAL temp = 0.0;
    INTEGER jc = 0;
    INTEGER ktemp = 0;
    INTEGER ctemp = 0;
    INTEGER nrp1 = 0;
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
    nlp1 = nl + 1;
    nlp2 = nl + 2;
    //
    if ((k < 1) || (k > n)) {
        info = -4;
    } else if (ldq < k) {
        info = -7;
    } else if (ldu < n) {
        info = -10;
    } else if (ldu2 < n) {
        info = -12;
    } else if (ldvt < m) {
        info = -14;
    } else if (ldvt2 < m) {
        info = -16;
    }
    if (info != 0) {
        Mxerbla("Rlasd3", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (k == 1) {
        d[1 - 1] = abs(z[1 - 1]);
        Rcopy(m, &vt2[(1 - 1)], ldvt2, &vt[(1 - 1)], ldvt);
        if (z[1 - 1] > zero) {
            Rcopy(n, &u2[(1 - 1)], 1, &u[(1 - 1)], 1);
        } else {
            for (i = 1; i <= n; i = i + 1) {
                u[(i - 1)] = -u2[(i - 1)];
            }
        }
        return;
    }
    //
    //     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
    //     be computed with high relative accuracy (barring over/underflow).
    //     This is a problem on machines without a guard digit in
    //     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
    //     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
    //     which on any of these machines zeros out the bottommost
    //     bit of DSIGMA(I) if it is 1; this makes the subsequent
    //     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
    //     occurs. On binary machines with a guard digit (almost all
    //     machines) it does not change DSIGMA(I) at all. On hexadecimal
    //     and decimal machines with a guard digit, it slightly
    //     changes the bottommost bits of DSIGMA(I). It does not account
    //     for hexadecimal or decimal machines without guard digits
    //     (we know of none). We use a subroutine call to compute
    //     2*DSIGMA(I) to prevent optimizing compilers from eliminating
    //     this code.
    //
    for (i = 1; i <= k; i = i + 1) {
        dsigma[i - 1] = Rlamc3(dsigma[i - 1], dsigma[i - 1]) - dsigma[i - 1];
    }
    //
    //     Keep a copy of Z.
    //
    Rcopy(k, z, 1, q, 1);
    //
    //     Normalize Z.
    //
    rho = Rnrm2(k, z, 1);
    Rlascl("G", 0, 0, rho, one, k, 1, z, k, info);
    rho = rho * rho;
    //
    //     Find the new singular values.
    //
    for (j = 1; j <= k; j = j + 1) {
        Rlasd4(k, j, dsigma, z, &u[(j - 1) * ldu], rho, d[j - 1], &vt[(j - 1) * ldvt], info);
        //
        //        If the zero finder fails, report the convergence failure.
        //
        if (info != 0) {
            return;
        }
    }
    //
    //     Compute updated Z.
    //
    for (i = 1; i <= k; i = i + 1) {
        z[i - 1] = u[(i - 1) + (k - 1) * ldu] * vt[(i - 1) + (k - 1) * ldvt];
        for (j = 1; j <= i - 1; j = j + 1) {
            z[i - 1] = z[i - 1] * (u[(i - 1) + (j - 1) * ldu] * vt[(i - 1) + (j - 1) * ldvt] / (dsigma[i - 1] - dsigma[j - 1]) / (dsigma[i - 1] + dsigma[j - 1]));
        }
        for (j = i; j <= k - 1; j = j + 1) {
            z[i - 1] = z[i - 1] * (u[(i - 1) + (j - 1) * ldu] * vt[(i - 1) + (j - 1) * ldvt] / (dsigma[i - 1] - dsigma[(j + 1) - 1]) / (dsigma[i - 1] + dsigma[(j + 1) - 1]));
        }
        z[i - 1] = sign(sqrt(abs(z[i - 1])), q[(i - 1)]);
    }
    //
    //     Compute left singular vectors of the modified diagonal matrix,
    //     and store related information for the right singular vectors.
    //
    for (i = 1; i <= k; i = i + 1) {
        vt[(i - 1) * ldvt] = z[1 - 1] / u[(i - 1) * ldu] / vt[(i - 1) * ldvt];
        u[(i - 1) * ldu] = negone;
        for (j = 2; j <= k; j = j + 1) {
            vt[(j - 1) + (i - 1) * ldvt] = z[j - 1] / u[(j - 1) + (i - 1) * ldu] / vt[(j - 1) + (i - 1) * ldvt];
            u[(j - 1) + (i - 1) * ldu] = dsigma[j - 1] * vt[(j - 1) + (i - 1) * ldvt];
        }
        temp = Rnrm2(k, &u[(i - 1) * ldu], (INTEGER)1);
        q[(i - 1) * ldq] = u[(i - 1) * ldu] / temp;
        for (j = 2; j <= k; j = j + 1) {
            jc = idxc[j - 1];
            q[(j - 1) + (i - 1) * ldq] = u[(jc - 1) + (i - 1) * ldu] / temp;
        }
    }
    //
    //     Update the left singular vector matrix.
    //
    if (k == 2) {
        Rgemm("N", "N", n, k, k, one, u2, ldu2, q, ldq, zero, u, ldu);
        goto statement_100;
    }
    if (ctot[1 - 1] > 0) {
        Rgemm("N", "N", nl, k, ctot[1 - 1], one, &u2[(2 - 1) * ldu2], ldu2, &q[(2 - 1)], ldq, zero, &u[(1 - 1)], ldu);
        if (ctot[3 - 1] > 0) {
            ktemp = 2 + ctot[1 - 1] + ctot[2 - 1];
            Rgemm("N", "N", nl, k, ctot[3 - 1], one, &u2[(ktemp - 1) * ldu2], ldu2, &q[(ktemp - 1)], ldq, one, &u[(1 - 1)], ldu);
        }
    } else if (ctot[3 - 1] > 0) {
        ktemp = 2 + ctot[1 - 1] + ctot[2 - 1];
        Rgemm("N", "N", nl, k, ctot[3 - 1], one, &u2[(ktemp - 1) * ldu2], ldu2, &q[(ktemp - 1)], ldq, zero, &u[(1 - 1)], ldu);
    } else {
        Rlacpy("F", nl, k, u2, ldu2, u, ldu);
    }
    Rcopy(k, &q[(1 - 1)], ldq, &u[(nlp1 - 1)], ldu);
    ktemp = 2 + ctot[1 - 1];
    ctemp = ctot[2 - 1] + ctot[3 - 1];
    Rgemm("N", "N", nr, k, ctemp, one, &u2[(nlp2 - 1) + (ktemp - 1) * ldu2], ldu2, &q[(ktemp - 1)], ldq, zero, &u[(nlp2 - 1)], ldu);
//
//     Generate the right singular vectors.
//
statement_100:
    for (i = 1; i <= k; i = i + 1) {
        temp = Rnrm2(k, &vt[(i - 1) * ldvt], 1);
        q[(i - 1)] = vt[(i - 1) * ldvt] / temp;
        for (j = 2; j <= k; j = j + 1) {
            jc = idxc[j - 1];
            q[(i - 1) + (j - 1) * ldq] = vt[(jc - 1) + (i - 1) * ldvt] / temp;
        }
    }
    //
    //     Update the right singular vector matrix.
    //
    if (k == 2) {
        Rgemm("N", "N", k, m, k, one, q, ldq, vt2, ldvt2, zero, vt, ldvt);
        return;
    }
    ktemp = 1 + ctot[1 - 1];
    Rgemm("N", "N", k, nlp1, ktemp, one, &q[(1 - 1)], ldq, &vt2[(1 - 1)], ldvt2, zero, &vt[(1 - 1)], ldvt);
    ktemp = 2 + ctot[1 - 1] + ctot[2 - 1];
    if (ktemp <= ldvt2) {
        Rgemm("N", "N", k, nlp1, ctot[3 - 1], one, &q[(ktemp - 1) * ldq], ldq, &vt2[(ktemp - 1)], ldvt2, one, &vt[(1 - 1)], ldvt);
    }
    //
    ktemp = ctot[1 - 1] + 1;
    nrp1 = nr + sqre;
    if (ktemp > 1) {
        for (i = 1; i <= k; i = i + 1) {
            q[(i - 1) + (ktemp - 1) * ldq] = q[(i - 1)];
        }
        for (i = nlp2; i <= m; i = i + 1) {
            vt2[(ktemp - 1) + (i - 1) * ldvt2] = vt2[(i - 1) * ldvt2];
        }
    }
    ctemp = 1 + ctot[2 - 1] + ctot[3 - 1];
    Rgemm("N", "N", k, nrp1, ctemp, one, &q[(ktemp - 1) * ldq], ldq, &vt2[(ktemp - 1) + (nlp2 - 1) * ldvt2], ldvt2, zero, &vt[(nlp2 - 1) * ldvt], ldvt);
    //
    //     End of Rlasd3
    //
}
