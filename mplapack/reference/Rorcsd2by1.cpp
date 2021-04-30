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

void Rorcsd2by1(const char *jobu1, const char *jobu2, const char *jobv1t, INTEGER const m, INTEGER const p, INTEGER const q, REAL *x11, INTEGER const ldx11, REAL *x21, INTEGER const ldx21, REAL *theta, REAL *u1, INTEGER const ldu1, REAL *u2, INTEGER const ldu2, REAL *v1t, INTEGER const ldv1t, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    //
    //  -- LAPACK computational routine (3.5.0) --
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Function ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test input arguments
    //
    info = 0;
    bool wantu1 = Mlsame(jobu1, "Y");
    bool wantu2 = Mlsame(jobu2, "Y");
    bool wantv1t = Mlsame(jobv1t, "Y");
    bool lquery = lwork == -1;
    //
    if (m < 0) {
        info = -4;
    } else if (p < 0 || p > m) {
        info = -5;
    } else if (q < 0 || q > m) {
        info = -6;
    } else if (ldx11 < max((INTEGER)1, p)) {
        info = -8;
    } else if (ldx21 < max((INTEGER)1, m - p)) {
        info = -10;
    } else if (wantu1 && ldu1 < max((INTEGER)1, p)) {
        info = -13;
    } else if (wantu2 && ldu2 < max((INTEGER)1, m - p)) {
        info = -15;
    } else if (wantv1t && ldv1t < max((INTEGER)1, q)) {
        info = -17;
    }
    //
    INTEGER r = min({p, m - p, q, m - q});
    //
    //     Compute workspace
    //
    //       WORK layout:
    //     |-------------------------------------------------------|
    //     | LWORKOPT (1)                                          |
    //     |-------------------------------------------------------|
    //     | PHI (MAX(1,R-1))                                      |
    //     |-------------------------------------------------------|
    //     | TAUP1 (MAX(1,P))                        | B11D (R)    |
    //     | TAUP2 (MAX(1,M-P))                      | B11E (R-1)  |
    //     | TAUQ1 (MAX(1,Q))                        | B12D (R)    |
    //     |-----------------------------------------| B12E (R-1)  |
    //     | Rorbdb WORK | Rorgqr WORK | Rorglq WORK | B21D (R)    |
    //     |             |             |             | B21E (R-1)  |
    //     |             |             |             | B22D (R)    |
    //     |             |             |             | B22E (R-1)  |
    //     |             |             |             | Rbbcsd WORK |
    //     |-------------------------------------------------------|
    //
    INTEGER iphi = 0;
    INTEGER ib11d = 0;
    INTEGER ib11e = 0;
    INTEGER ib12d = 0;
    INTEGER ib12e = 0;
    INTEGER ib21d = 0;
    INTEGER ib21e = 0;
    INTEGER ib22d = 0;
    INTEGER ib22e = 0;
    INTEGER ibbcsd = 0;
    INTEGER itaup1 = 0;
    INTEGER itaup2 = 0;
    INTEGER itauq1 = 0;
    INTEGER iorbdb = 0;
    INTEGER iorgqr = 0;
    INTEGER iorglq = 0;
    INTEGER lorgqrmin = 0;
    INTEGER lorgqropt = 0;
    INTEGER lorglqmin = 0;
    INTEGER lorglqopt = 0;
    REAL dum1[1];
    INTEGER childinfo = 0;
    INTEGER lorbdb = 0;
    REAL dum2[1];
    INTEGER lbbcsd = 0;
    INTEGER lworkmin = 0;
    INTEGER lworkopt = 0;
    if (info == 0) {
        iphi = 2;
        ib11d = iphi + max((INTEGER)1, r - 1);
        ib11e = ib11d + max((INTEGER)1, r);
        ib12d = ib11e + max((INTEGER)1, r - 1);
        ib12e = ib12d + max((INTEGER)1, r);
        ib21d = ib12e + max((INTEGER)1, r - 1);
        ib21e = ib21d + max((INTEGER)1, r);
        ib22d = ib21e + max((INTEGER)1, r - 1);
        ib22e = ib22d + max((INTEGER)1, r);
        ibbcsd = ib22e + max((INTEGER)1, r - 1);
        itaup1 = iphi + max((INTEGER)1, r - 1);
        itaup2 = itaup1 + max((INTEGER)1, p);
        itauq1 = itaup2 + max((INTEGER)1, m - p);
        iorbdb = itauq1 + max((INTEGER)1, q);
        iorgqr = itauq1 + max((INTEGER)1, q);
        iorglq = itauq1 + max((INTEGER)1, q);
        lorgqrmin = 1;
        lorgqropt = 1;
        lorglqmin = 1;
        lorglqopt = 1;
        if (r == q) {
            Rorbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, dum1, dum1, dum1, dum1, work, -1, childinfo);
            lorbdb = castINTEGER(work[1 - 1]);
            if (wantu1 && p > 0) {
                Rorgqr(p, p, q, u1, ldu1, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, p);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantu2 && m - p > 0) {
                Rorgqr(m - p, m - p, q, u2, ldu2, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, m - p);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantv1t && q > 0) {
                Rorglq(q - 1, q - 1, q - 1, v1t, ldv1t, dum1, &work[1 - 1], -1, childinfo);
                lorglqmin = max(lorglqmin, q - 1);
                lorglqopt = max(lorglqopt, castINTEGER(work[1 - 1]));
            }
            Rbbcsd(jobu1, jobu2, jobv1t, "N", "N", m, p, q, theta, dum1, u1, ldu1, u2, ldu2, v1t, ldv1t, dum2, 1, dum1, dum1, dum1, dum1, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lbbcsd = castINTEGER(work[1 - 1]);
        } else if (r == p) {
            Rorbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lorbdb = castINTEGER(work[1 - 1]);
            if (wantu1 && p > 0) {
                Rorgqr(p - 1, p - 1, p - 1, &u1[(2 - 1) + (2 - 1) * ldu1], ldu1, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, p - 1);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantu2 && m - p > 0) {
                Rorgqr(m - p, m - p, q, u2, ldu2, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, m - p);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantv1t && q > 0) {
                Rorglq(q, q, r, v1t, ldv1t, dum1, &work[1 - 1], -1, childinfo);
                lorglqmin = max(lorglqmin, q);
                lorglqopt = max(lorglqopt, castINTEGER(work[1 - 1]));
            }
            Rbbcsd(jobv1t, "N", jobu1, jobu2, "T", m, q, p, theta, dum1, v1t, ldv1t, dum2, 1, u1, ldu1, u2, ldu2, dum1, dum1, dum1, dum1, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lbbcsd = castINTEGER(work[1 - 1]);
        } else if (r == m - p) {
            Rorbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lorbdb = castINTEGER(work[1 - 1]);
            if (wantu1 && p > 0) {
                Rorgqr(p, p, q, u1, ldu1, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, p);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantu2 && m - p > 0) {
                Rorgqr(m - p - 1, m - p - 1, m - p - 1, &u2[(2 - 1) + (2 - 1) * ldu2], ldu2, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, m - p - 1);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantv1t && q > 0) {
                Rorglq(q, q, r, v1t, ldv1t, dum1, &work[1 - 1], -1, childinfo);
                lorglqmin = max(lorglqmin, q);
                lorglqopt = max(lorglqopt, castINTEGER(work[1 - 1]));
            }
            Rbbcsd("N", jobv1t, jobu2, jobu1, "T", m, m - q, m - p, theta, dum1, dum2, 1, v1t, ldv1t, u2, ldu2, u1, ldu1, dum1, dum1, dum1, dum1, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lbbcsd = castINTEGER(work[1 - 1]);
        } else {
            Rorbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, dum1, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lorbdb = m + castINTEGER(work[1 - 1]);
            if (wantu1 && p > 0) {
                Rorgqr(p, p, m - q, u1, ldu1, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, p);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantu2 && m - p > 0) {
                Rorgqr(m - p, m - p, m - q, u2, ldu2, dum1, &work[1 - 1], -1, childinfo);
                lorgqrmin = max(lorgqrmin, m - p);
                lorgqropt = max(lorgqropt, castINTEGER(work[1 - 1]));
            }
            if (wantv1t && q > 0) {
                Rorglq(q, q, q, v1t, ldv1t, dum1, &work[1 - 1], -1, childinfo);
                lorglqmin = max(lorglqmin, q);
                lorglqopt = max(lorglqopt, castINTEGER(work[1 - 1]));
            }
            Rbbcsd(jobu2, jobu1, "N", jobv1t, "N", m, m - p, m - q, theta, dum1, u2, ldu2, u1, ldu1, dum2, 1, v1t, ldv1t, dum1, dum1, dum1, dum1, dum1, dum1, dum1, dum1, &work[1 - 1], -1, childinfo);
            lbbcsd = castINTEGER(work[1 - 1]);
        }
        lworkmin = max({iorbdb + lorbdb - 1, iorgqr + lorgqrmin - 1, iorglq + lorglqmin - 1, ibbcsd + lbbcsd - 1});
        lworkopt = max({iorbdb + lorbdb - 1, iorgqr + lorgqropt - 1, iorglq + lorglqopt - 1, ibbcsd + lbbcsd - 1});
        work[1 - 1] = lworkopt;
        if (lwork < lworkmin && !lquery) {
            info = -19;
        }
    }
    if (info != 0) {
        Mxerbla("Rorcsd2by1", -info);
        return;
    } else if (lquery) {
        return;
    }
    INTEGER lorgqr = lwork - iorgqr + 1;
    INTEGER lorglq = lwork - iorglq + 1;
    //
    //     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
    //     in which R = MIN(P,M-P,Q,M-Q)
    //
    const REAL one = 1.0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    INTEGER i = 0;
    if (r == q) {
        //
        //        Case 1: R = Q
        //
        //        Simultaneously bidiagonalize X11 and X21
        //
        Rorbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, &work[iphi - 1], &work[itaup1 - 1], &work[itaup2 - 1], &work[itauq1 - 1], &work[iorbdb - 1], lorbdb, childinfo);
        //
        //        Accumulate Householder reflectors
        //
        if (wantu1 && p > 0) {
            Rlacpy("L", p, q, x11, ldx11, u1, ldu1);
            Rorgqr(p, p, q, u1, ldu1, &work[itaup1 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantu2 && m - p > 0) {
            Rlacpy("L", m - p, q, x21, ldx21, u2, ldu2);
            Rorgqr(m - p, m - p, q, u2, ldu2, &work[itaup2 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantv1t && q > 0) {
            v1t[(1 - 1)] = one;
            for (j = 2; j <= q; j = j + 1) {
                v1t[(j - 1) * ldv1t] = zero;
                v1t[(j - 1)] = zero;
            }
            Rlacpy("U", q - 1, q - 1, &x21[(2 - 1) * ldx21], ldx21, &v1t[(2 - 1) + (2 - 1) * ldv1t], ldv1t);
            Rorglq(q - 1, q - 1, q - 1, &v1t[(2 - 1) + (2 - 1) * ldv1t], ldv1t, &work[itauq1 - 1], &work[iorglq - 1], lorglq, childinfo);
        }
        //
        //        Simultaneously diagonalize X11 and X21.
        //
        Rbbcsd(jobu1, jobu2, jobv1t, "N", "N", m, p, q, theta, &work[iphi - 1], u1, ldu1, u2, ldu2, v1t, ldv1t, dum2, 1, &work[ib11d - 1], &work[ib11e - 1], &work[ib12d - 1], &work[ib12e - 1], &work[ib21d - 1], &work[ib21e - 1], &work[ib22d - 1], &work[ib22e - 1], &work[ibbcsd - 1], lbbcsd, childinfo);
        //
        //        Permute rows and columns to place zero submatrices in
        //        preferred positions
        //
        if (q > 0 && wantu2) {
            for (i = 1; i <= q; i = i + 1) {
                iwork[i - 1] = m - p - q + i;
            }
            for (i = q + 1; i <= m - p; i = i + 1) {
                iwork[i - 1] = i - q;
            }
            Rlapmt(false, m - p, m - p, u2, ldu2, iwork);
        }
    } else if (r == p) {
        //
        //        Case 2: R = P
        //
        //        Simultaneously bidiagonalize X11 and X21
        //
        Rorbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, &work[iphi - 1], &work[itaup1 - 1], &work[itaup2 - 1], &work[itauq1 - 1], &work[iorbdb - 1], lorbdb, childinfo);
        //
        //        Accumulate Householder reflectors
        //
        if (wantu1 && p > 0) {
            u1[(1 - 1)] = one;
            for (j = 2; j <= p; j = j + 1) {
                u1[(j - 1) * ldu1] = zero;
                u1[(j - 1)] = zero;
            }
            Rlacpy("L", p - 1, p - 1, &x11[(2 - 1)], ldx11, &u1[(2 - 1) + (2 - 1) * ldu1], ldu1);
            Rorgqr(p - 1, p - 1, p - 1, &u1[(2 - 1) + (2 - 1) * ldu1], ldu1, &work[itaup1 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantu2 && m - p > 0) {
            Rlacpy("L", m - p, q, x21, ldx21, u2, ldu2);
            Rorgqr(m - p, m - p, q, u2, ldu2, &work[itaup2 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantv1t && q > 0) {
            Rlacpy("U", p, q, x11, ldx11, v1t, ldv1t);
            Rorglq(q, q, r, v1t, ldv1t, &work[itauq1 - 1], &work[iorglq - 1], lorglq, childinfo);
        }
        //
        //        Simultaneously diagonalize X11 and X21.
        //
        Rbbcsd(jobv1t, "N", jobu1, jobu2, "T", m, q, p, theta, &work[iphi - 1], v1t, ldv1t, dum2, 1, u1, ldu1, u2, ldu2, &work[ib11d - 1], &work[ib11e - 1], &work[ib12d - 1], &work[ib12e - 1], &work[ib21d - 1], &work[ib21e - 1], &work[ib22d - 1], &work[ib22e - 1], &work[ibbcsd - 1], lbbcsd, childinfo);
        //
        //        Permute rows and columns to place identity submatrices in
        //        preferred positions
        //
        if (q > 0 && wantu2) {
            for (i = 1; i <= q; i = i + 1) {
                iwork[i - 1] = m - p - q + i;
            }
            for (i = q + 1; i <= m - p; i = i + 1) {
                iwork[i - 1] = i - q;
            }
            Rlapmt(false, m - p, m - p, u2, ldu2, iwork);
        }
    } else if (r == m - p) {
        //
        //        Case 3: R = M-P
        //
        //        Simultaneously bidiagonalize X11 and X21
        //
        Rorbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, &work[iphi - 1], &work[itaup1 - 1], &work[itaup2 - 1], &work[itauq1 - 1], &work[iorbdb - 1], lorbdb, childinfo);
        //
        //        Accumulate Householder reflectors
        //
        if (wantu1 && p > 0) {
            Rlacpy("L", p, q, x11, ldx11, u1, ldu1);
            Rorgqr(p, p, q, u1, ldu1, &work[itaup1 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantu2 && m - p > 0) {
            u2[(1 - 1)] = one;
            for (j = 2; j <= m - p; j = j + 1) {
                u2[(j - 1) * ldu2] = zero;
                u2[(j - 1)] = zero;
            }
            Rlacpy("L", m - p - 1, m - p - 1, &x21[(2 - 1)], ldx21, &u2[(2 - 1) + (2 - 1) * ldu2], ldu2);
            Rorgqr(m - p - 1, m - p - 1, m - p - 1, &u2[(2 - 1) + (2 - 1) * ldu2], ldu2, &work[itaup2 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantv1t && q > 0) {
            Rlacpy("U", m - p, q, x21, ldx21, v1t, ldv1t);
            Rorglq(q, q, r, v1t, ldv1t, &work[itauq1 - 1], &work[iorglq - 1], lorglq, childinfo);
        }
        //
        //        Simultaneously diagonalize X11 and X21.
        //
        Rbbcsd("N", jobv1t, jobu2, jobu1, "T", m, m - q, m - p, theta, &work[iphi - 1], dum2, 1, v1t, ldv1t, u2, ldu2, u1, ldu1, &work[ib11d - 1], &work[ib11e - 1], &work[ib12d - 1], &work[ib12e - 1], &work[ib21d - 1], &work[ib21e - 1], &work[ib22d - 1], &work[ib22e - 1], &work[ibbcsd - 1], lbbcsd, childinfo);
        //
        //        Permute rows and columns to place identity submatrices in
        //        preferred positions
        //
        if (q > r) {
            for (i = 1; i <= r; i = i + 1) {
                iwork[i - 1] = q - r + i;
            }
            for (i = r + 1; i <= q; i = i + 1) {
                iwork[i - 1] = i - r;
            }
            if (wantu1) {
                Rlapmt(false, p, q, u1, ldu1, iwork);
            }
            if (wantv1t) {
                Rlapmr(false, q, q, v1t, ldv1t, iwork);
            }
        }
    } else {
        //
        //        Case 4: R = M-Q
        //
        //        Simultaneously bidiagonalize X11 and X21
        //
        Rorbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, &work[iphi - 1], &work[itaup1 - 1], &work[itaup2 - 1], &work[itauq1 - 1], &work[iorbdb - 1], &work[(iorbdb + m) - 1], lorbdb - m, childinfo);
        //
        //        Accumulate Householder reflectors
        //
        if (wantu2 && m - p > 0) {
            Rcopy(m - p, &work[(iorbdb + p) - 1], 1, u2, 1);
        }
        if (wantu1 && p > 0) {
            Rcopy(p, &work[iorbdb - 1], 1, u1, 1);
            for (j = 2; j <= p; j = j + 1) {
                u1[(j - 1) * ldu1] = zero;
            }
            Rlacpy("L", p - 1, m - q - 1, &x11[(2 - 1)], ldx11, &u1[(2 - 1) + (2 - 1) * ldu1], ldu1);
            Rorgqr(p, p, m - q, u1, ldu1, &work[itaup1 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantu2 && m - p > 0) {
            for (j = 2; j <= m - p; j = j + 1) {
                u2[(j - 1) * ldu2] = zero;
            }
            Rlacpy("L", m - p - 1, m - q - 1, &x21[(2 - 1)], ldx21, &u2[(2 - 1) + (2 - 1) * ldu2], ldu2);
            Rorgqr(m - p, m - p, m - q, u2, ldu2, &work[itaup2 - 1], &work[iorgqr - 1], lorgqr, childinfo);
        }
        if (wantv1t && q > 0) {
            Rlacpy("U", m - q, q, x21, ldx21, v1t, ldv1t);
            Rlacpy("U", p - (m - q), q - (m - q), &x11[((m - q + 1) - 1) + ((m - q + 1) - 1) * ldx11], ldx11, &v1t[((m - q + 1) - 1) + ((m - q + 1) - 1) * ldv1t], ldv1t);
            Rlacpy("U", -p + q, q - p, &x21[((m - q + 1) - 1) + ((p + 1) - 1) * ldx21], ldx21, &v1t[((p + 1) - 1) + ((p + 1) - 1) * ldv1t], ldv1t);
            Rorglq(q, q, q, v1t, ldv1t, &work[itauq1 - 1], &work[iorglq - 1], lorglq, childinfo);
        }
        //
        //        Simultaneously diagonalize X11 and X21.
        //
        Rbbcsd(jobu2, jobu1, "N", jobv1t, "N", m, m - p, m - q, theta, &work[iphi - 1], u2, ldu2, u1, ldu1, dum2, 1, v1t, ldv1t, &work[ib11d - 1], &work[ib11e - 1], &work[ib12d - 1], &work[ib12e - 1], &work[ib21d - 1], &work[ib21e - 1], &work[ib22d - 1], &work[ib22e - 1], &work[ibbcsd - 1], lbbcsd, childinfo);
        //
        //        Permute rows and columns to place identity submatrices in
        //        preferred positions
        //
        if (p > r) {
            for (i = 1; i <= r; i = i + 1) {
                iwork[i - 1] = p - r + i;
            }
            for (i = r + 1; i <= p; i = i + 1) {
                iwork[i - 1] = i - r;
            }
            if (wantu1) {
                Rlapmt(false, p, p, u1, ldu1, iwork);
            }
            if (wantv1t) {
                Rlapmr(false, p, q, v1t, ldv1t, iwork);
            }
        }
    }
    //
    //     End of Rorcsd2by1
    //
}
