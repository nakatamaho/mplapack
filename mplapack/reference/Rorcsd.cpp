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

void Rorcsd(const char *jobu1, const char *jobu2, const char *jobv1t, const char *jobv2t, const char *trans, const char *signs, INTEGER const &m, INTEGER const &p, INTEGER const &q, REAL *x11, INTEGER const &ldx11, REAL *x12, INTEGER const &ldx12, REAL *x21, INTEGER const &ldx21, REAL *x22, INTEGER const &ldx22, REAL *theta, REAL *u1, INTEGER const &ldu1, REAL *u2, INTEGER const &ldu2, REAL *v1t, INTEGER const &ldv1t, REAL *v2t, INTEGER const &ldv2t, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER &info) {
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
    //  ===================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions
    //     ..
    //     .. Executable Statements ..
    //
    //     Test input arguments
    //
    info = 0;
    bool wantu1 = Mlsame(jobu1, "Y");
    bool wantu2 = Mlsame(jobu2, "Y");
    bool wantv1t = Mlsame(jobv1t, "Y");
    bool wantv2t = Mlsame(jobv2t, "Y");
    bool colmajor = !Mlsame(trans, "T");
    bool defaultsigns = !Mlsame(signs, "O");
    bool lquery = lwork == -1;
    if (m < 0) {
        info = -7;
    } else if (p < 0 || p > m) {
        info = -8;
    } else if (q < 0 || q > m) {
        info = -9;
    } else if (colmajor && ldx11 < max((INTEGER)1, p)) {
        info = -11;
    } else if (!colmajor && ldx11 < max((INTEGER)1, q)) {
        info = -11;
    } else if (colmajor && ldx12 < max((INTEGER)1, p)) {
        info = -13;
    } else if (!colmajor && ldx12 < max((INTEGER)1, m - q)) {
        info = -13;
    } else if (colmajor && ldx21 < max((INTEGER)1, m - p)) {
        info = -15;
    } else if (!colmajor && ldx21 < max((INTEGER)1, q)) {
        info = -15;
    } else if (colmajor && ldx22 < max((INTEGER)1, m - p)) {
        info = -17;
    } else if (!colmajor && ldx22 < max((INTEGER)1, m - q)) {
        info = -17;
    } else if (wantu1 && ldu1 < p) {
        info = -20;
    } else if (wantu2 && ldu2 < m - p) {
        info = -22;
    } else if (wantv1t && ldv1t < q) {
        info = -24;
    } else if (wantv2t && ldv2t < m - q) {
        info = -26;
    }
    //
    //     Work with transpose if convenient
    //
    str<1> transt = char0;
    str<1> signst = char0;
    if (info == 0 && min(p, m - p) < min(q, m - q)) {
        if (colmajor) {
            transt = "T";
        } else {
            transt = "N";
        }
        if (defaultsigns) {
            signst = "O";
        } else {
            signst = "D";
        }
        Rorcsd(jobv1t, jobv2t, jobu1, jobu2, transt, signst, m, q, p, x11, ldx11, x21, ldx21, x12, ldx12, x22, ldx22, theta, v1t, ldv1t, v2t, ldv2t, u1, ldu1, u2, ldu2, work, lwork, iwork, info);
        return;
    }
    //
    //     Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
    //     convenient
    //
    if (info == 0 && m - q < q) {
        if (defaultsigns) {
            signst = "O";
        } else {
            signst = "D";
        }
        Rorcsd(jobu2, jobu1, jobv2t, jobv1t, trans, signst, m, m - p, m - q, x22, ldx22, x21, ldx21, x12, ldx12, x11, ldx11, theta, u2, ldu2, u1, ldu1, v2t, ldv2t, v1t, ldv1t, work, lwork, iwork, info);
        return;
    }
    //
    //     Compute workspace
    //
    INTEGER iphi = 0;
    INTEGER itaup1 = 0;
    INTEGER itaup2 = 0;
    INTEGER itauq1 = 0;
    INTEGER itauq2 = 0;
    INTEGER iorgqr = 0;
    INTEGER childinfo = 0;
    INTEGER lorgqrworkopt = 0;
    INTEGER lorgqrworkmin = 0;
    INTEGER iorglq = 0;
    INTEGER lorglqworkopt = 0;
    INTEGER lorglqworkmin = 0;
    INTEGER iorbdb = 0;
    INTEGER lorbdbworkopt = 0;
    INTEGER lorbdbworkmin = 0;
    INTEGER ib11d = 0;
    INTEGER ib11e = 0;
    INTEGER ib12d = 0;
    INTEGER ib12e = 0;
    INTEGER ib21d = 0;
    INTEGER ib21e = 0;
    INTEGER ib22d = 0;
    INTEGER ib22e = 0;
    INTEGER ibbcsd = 0;
    INTEGER lbbcsdworkopt = 0;
    INTEGER lbbcsdworkmin = 0;
    INTEGER lworkopt = 0;
    INTEGER lworkmin = 0;
    INTEGER lorgqrwork = 0;
    INTEGER lorglqwork = 0;
    INTEGER lorbdbwork = 0;
    INTEGER lbbcsdwork = 0;
    if (info == 0) {
        //
        iphi = 2;
        itaup1 = iphi + max((INTEGER)1, q - 1);
        itaup2 = itaup1 + max((INTEGER)1, p);
        itauq1 = itaup2 + max((INTEGER)1, m - p);
        itauq2 = itauq1 + max((INTEGER)1, q);
        iorgqr = itauq2 + max((INTEGER)1, m - q);
        Rorgqr(m - q, m - q, m - q, u1, max((INTEGER)1, m - q), u1, work, -1, childinfo);
        lorgqrworkopt = INTEGER(work[1 - 1]);
        lorgqrworkmin = max((INTEGER)1, m - q);
        iorglq = itauq2 + max((INTEGER)1, m - q);
        Rorglq(m - q, m - q, m - q, u1, max((INTEGER)1, m - q), u1, work, -1, childinfo);
        lorglqworkopt = INTEGER(work[1 - 1]);
        lorglqworkmin = max((INTEGER)1, m - q);
        iorbdb = itauq2 + max((INTEGER)1, m - q);
        Rorbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, v1t, u1, u2, v1t, v2t, work, -1, childinfo);
        lorbdbworkopt = INTEGER(work[1 - 1]);
        lorbdbworkmin = lorbdbworkopt;
        ib11d = itauq2 + max((INTEGER)1, m - q);
        ib11e = ib11d + max((INTEGER)1, q);
        ib12d = ib11e + max((INTEGER)1, q - 1);
        ib12e = ib12d + max((INTEGER)1, q);
        ib21d = ib12e + max((INTEGER)1, q - 1);
        ib21e = ib21d + max((INTEGER)1, q);
        ib22d = ib21e + max((INTEGER)1, q - 1);
        ib22e = ib22d + max((INTEGER)1, q);
        ibbcsd = ib22e + max((INTEGER)1, q - 1);
        Rbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, u1, u1, u1, u1, u1, u1, u1, u1, work, -1, childinfo);
        lbbcsdworkopt = INTEGER(work[1 - 1]);
        lbbcsdworkmin = lbbcsdworkopt;
        lworkopt = max(iorgqr + lorgqrworkopt, iorglq + lorglqworkopt, iorbdb + lorbdbworkopt, ibbcsd + lbbcsdworkopt) - 1;
        lworkmin = max(iorgqr + lorgqrworkmin, iorglq + lorglqworkmin, iorbdb + lorbdbworkopt, ibbcsd + lbbcsdworkmin) - 1;
        work[1 - 1] = max(lworkopt, lworkmin);
        //
        if (lwork < lworkmin && !lquery) {
            info = -22;
        } else {
            lorgqrwork = lwork - iorgqr + 1;
            lorglqwork = lwork - iorglq + 1;
            lorbdbwork = lwork - iorbdb + 1;
            lbbcsdwork = lwork - ibbcsd + 1;
        }
    }
    //
    //     Abort if any illegal arguments
    //
    if (info != 0) {
        Mxerbla("Rorcsd", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Transform to bidiagonal block form
    //
    Rorbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, work[iphi - 1], work[itaup1 - 1], work[itaup2 - 1], work[itauq1 - 1], work[itauq2 - 1], work[iorbdb - 1], lorbdbwork, childinfo);
    //
    //     Accumulate Householder reflectors
    //
    const REAL one = 1.0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    if (colmajor) {
        if (wantu1 && p > 0) {
            Rlacpy("L", p, q, x11, ldx11, u1, ldu1);
            Rorgqr(p, p, q, u1, ldu1, work[itaup1 - 1], work[iorgqr - 1], lorgqrwork, info);
        }
        if (wantu2 && m - p > 0) {
            Rlacpy("L", m - p, q, x21, ldx21, u2, ldu2);
            Rorgqr(m - p, m - p, q, u2, ldu2, work[itaup2 - 1], work[iorgqr - 1], lorgqrwork, info);
        }
        if (wantv1t && q > 0) {
            Rlacpy("U", q - 1, q - 1, x11[(2 - 1) * ldx11], ldx11, v1t[(2 - 1) + (2 - 1) * ldv1t], ldv1t);
            v1t[(1 - 1)] = one;
            for (j = 2; j <= q; j = j + 1) {
                v1t[(j - 1) * ldv1t] = zero;
                v1t[(j - 1)] = zero;
            }
            Rorglq(q - 1, q - 1, q - 1, v1t[(2 - 1) + (2 - 1) * ldv1t], ldv1t, work[itauq1 - 1], work[iorglq - 1], lorglqwork, info);
        }
        if (wantv2t && m - q > 0) {
            Rlacpy("U", p, m - q, x12, ldx12, v2t, ldv2t);
            if (m - p > q) {
                Rlacpy("U", m - p - q, m - p - q, x22[((q + 1) - 1) + ((p + 1) - 1) * ldx22], ldx22, v2t[((p + 1) - 1) + ((p + 1) - 1) * ldv2t], ldv2t);
            }
            if (m > q) {
                Rorglq(m - q, m - q, m - q, v2t, ldv2t, work[itauq2 - 1], work[iorglq - 1], lorglqwork, info);
            }
        }
    } else {
        if (wantu1 && p > 0) {
            Rlacpy("U", q, p, x11, ldx11, u1, ldu1);
            Rorglq(p, p, q, u1, ldu1, work[itaup1 - 1], work[iorglq - 1], lorglqwork, info);
        }
        if (wantu2 && m - p > 0) {
            Rlacpy("U", q, m - p, x21, ldx21, u2, ldu2);
            Rorglq(m - p, m - p, q, u2, ldu2, work[itaup2 - 1], work[iorglq - 1], lorglqwork, info);
        }
        if (wantv1t && q > 0) {
            Rlacpy("L", q - 1, q - 1, x11[(2 - 1)], ldx11, v1t[(2 - 1) + (2 - 1) * ldv1t], ldv1t);
            v1t[(1 - 1)] = one;
            for (j = 2; j <= q; j = j + 1) {
                v1t[(j - 1) * ldv1t] = zero;
                v1t[(j - 1)] = zero;
            }
            Rorgqr(q - 1, q - 1, q - 1, v1t[(2 - 1) + (2 - 1) * ldv1t], ldv1t, work[itauq1 - 1], work[iorgqr - 1], lorgqrwork, info);
        }
        if (wantv2t && m - q > 0) {
            Rlacpy("L", m - q, p, x12, ldx12, v2t, ldv2t);
            Rlacpy("L", m - p - q, m - p - q, x22[((p + 1) - 1) + ((q + 1) - 1) * ldx22], ldx22, v2t[((p + 1) - 1) + ((p + 1) - 1) * ldv2t], ldv2t);
            Rorgqr(m - q, m - q, m - q, v2t, ldv2t, work[itauq2 - 1], work[iorgqr - 1], lorgqrwork, info);
        }
    }
    //
    //     Compute the CSD of the matrix in bidiagonal-block form
    //
    Rbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, work[iphi - 1], u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work[ib11d - 1], work[ib11e - 1], work[ib12d - 1], work[ib12e - 1], work[ib21d - 1], work[ib21e - 1], work[ib22d - 1], work[ib22e - 1], work[ibbcsd - 1], lbbcsdwork, info);
    //
    //     Permute rows and columns to place identity submatrices in top-
    //     left corner of (1,1)-block and/or bottom-right corner of (1,2)-
    //     block and/or bottom-right corner of (2,1)-block and/or top-left
    //     corner of (2,2)-block
    //
    INTEGER i = 0;
    if (q > 0 && wantu2) {
        for (i = 1; i <= q; i = i + 1) {
            iwork[i - 1] = m - p - q + i;
        }
        for (i = q + 1; i <= m - p; i = i + 1) {
            iwork[i - 1] = i - q;
        }
        if (colmajor) {
            Rlapmt(false, m - p, m - p, u2, ldu2, iwork);
        } else {
            Rlapmr(false, m - p, m - p, u2, ldu2, iwork);
        }
    }
    if (m > 0 && wantv2t) {
        for (i = 1; i <= p; i = i + 1) {
            iwork[i - 1] = m - p - q + i;
        }
        for (i = p + 1; i <= m - q; i = i + 1) {
            iwork[i - 1] = i - p;
        }
        if (!colmajor) {
            Rlapmt(false, m - q, m - q, v2t, ldv2t, iwork);
        } else {
            Rlapmr(false, m - q, m - q, v2t, ldv2t, iwork);
        }
    }
    //
    //     End Rorcsd
    //
}
