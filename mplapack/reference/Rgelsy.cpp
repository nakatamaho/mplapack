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

void Rgelsy(INTEGER const m, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER *jpvt, REAL const rcond, INTEGER &rank, REAL *work, INTEGER const lwork, INTEGER &info) {
    INTEGER mn = 0;
    INTEGER ismin = 0;
    INTEGER ismax = 0;
    bool lquery = false;
    INTEGER lwkmin = 0;
    INTEGER lwkopt = 0;
    INTEGER nb1 = 0;
    INTEGER nb2 = 0;
    INTEGER nb3 = 0;
    INTEGER nb4 = 0;
    INTEGER nb = 0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    INTEGER iascl = 0;
    const REAL zero = 0.0;
    REAL bnrm = 0.0;
    INTEGER ibscl = 0;
    REAL wsize = 0.0;
    REAL smax = 0.0;
    REAL smin = 0.0;
    INTEGER i = 0;
    const INTEGER imin = 2;
    REAL sminpr = 0.0;
    REAL s1 = 0.0;
    REAL c1 = 0.0;
    const INTEGER imax = 1;
    REAL smaxpr = 0.0;
    REAL s2 = 0.0;
    REAL c2 = 0.0;
    INTEGER j = 0;
    //
    //  -- LAPACK driver routine --
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
    mn = min(m, n);
    ismin = mn + 1;
    ismax = 2 * mn + 1;
    //
    //     Test the input arguments.
    //
    info = 0;
    lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldb < max({(INTEGER)1, m, n})) {
        info = -7;
    }
    //
    //     Figure out optimal block size
    //
    if (info == 0) {
        if (mn == 0 || nrhs == 0) {
            lwkmin = 1;
            lwkopt = 1;
        } else {
            nb1 = iMlaenv(1, "Rgeqrf", " ", m, n, -1, -1);
            nb2 = iMlaenv(1, "RgerQF", " ", m, n, -1, -1);
            nb3 = iMlaenv(1, "Rormqr", " ", m, n, nrhs, -1);
            nb4 = iMlaenv(1, "Rormrq", " ", m, n, nrhs, -1);
            nb = max({nb1, nb2, nb3, nb4});
            lwkmin = mn + max({2 * mn, n + 1, mn + nrhs});
            lwkopt = max({lwkmin, mn + 2 * n + nb * (n + 1), 2 * mn + nb * nrhs});
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < lwkmin && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgelsy", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (mn == 0 || nrhs == 0) {
        rank = 0;
        return;
    }
    //
    //     Get machine parameters
    //
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Scale A, B if max entries outside range [SMLNUM,BIGNUM]
    //
    anrm = Rlange("M", m, n, a, lda, work);
    iascl = 0;
    if (anrm > zero && anrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Rlascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        iascl = 1;
    } else if (anrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Rlascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
        iascl = 2;
    } else if (anrm == zero) {
        //
        //        Matrix all zero. Return zero solution.
        //
        Rlaset("F", max(m, n), nrhs, zero, zero, b, ldb);
        rank = 0;
        goto statement_70;
    }
    //
    bnrm = Rlange("M", m, nrhs, b, ldb, work);
    ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Rlascl("G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Rlascl("G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
        ibscl = 2;
    }
    //
    //     Compute QR factorization with column pivoting of A:
    //        A * P = Q * R
    //
    Rgeqp3(m, n, a, lda, jpvt, &work[1 - 1], &work[(mn + 1) - 1], lwork - mn, info);
    wsize = mn + work[(mn + 1) - 1];
    //
    //     workspace: MN+2*N+NB*(N+1).
    //     Details of Householder rotations stored in WORK(1:MN).
    //
    //     Determine RANK using incremental condition estimation
    //
    work[ismin - 1] = one;
    work[ismax - 1] = one;
    smax = abs(a[(1 - 1)]);
    smin = smax;
    if (abs(a[(1 - 1)]) == zero) {
        rank = 0;
        Rlaset("F", max(m, n), nrhs, zero, zero, b, ldb);
        goto statement_70;
    } else {
        rank = 1;
    }
//
statement_10:
    if (rank < mn) {
        i = rank + 1;
        Rlaic1(imin, rank, &work[ismin - 1], smin, &a[(i - 1) * lda], a[(i - 1) + (i - 1) * lda], sminpr, s1, c1);
        Rlaic1(imax, rank, &work[ismax - 1], smax, &a[(i - 1) * lda], a[(i - 1) + (i - 1) * lda], smaxpr, s2, c2);
        //
        if (smaxpr * rcond <= sminpr) {
            for (i = 1; i <= rank; i = i + 1) {
                work[(ismin + i - 1) - 1] = s1 * work[(ismin + i - 1) - 1];
                work[(ismax + i - 1) - 1] = s2 * work[(ismax + i - 1) - 1];
            }
            work[(ismin + rank) - 1] = c1;
            work[(ismax + rank) - 1] = c2;
            smin = sminpr;
            smax = smaxpr;
            rank++;
            goto statement_10;
        }
    }
    //
    //     workspace: 3*MN.
    //
    //     Logically partition R = [ R11 R12 ]
    //                             [  0  R22 ]
    //     where R11 = R(1:RANK,1:RANK)
    //
    //     [R11,R12] = [ T11, 0 ] * Y
    //
    if (rank < n) {
        Rtzrzf(rank, n, a, lda, &work[(mn + 1) - 1], &work[(2 * mn + 1) - 1], lwork - 2 * mn, info);
    }
    //
    //     workspace: 2*MN.
    //     Details of Householder rotations stored in WORK(MN+1:2*MN)
    //
    //     B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
    //
    Rormqr("Left", "Transpose", m, nrhs, mn, a, lda, &work[1 - 1], b, ldb, &work[(2 * mn + 1) - 1], lwork - 2 * mn, info);
    wsize = max(wsize, REAL(2 * mn + work[(2 * mn + 1) - 1]));
    //
    //     workspace: 2*MN+NB*NRHS.
    //
    //     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
    //
    Rtrsm("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, one, a, lda, b, ldb);
    //
    for (j = 1; j <= nrhs; j = j + 1) {
        for (i = rank + 1; i <= n; i = i + 1) {
            b[(i - 1) + (j - 1) * ldb] = zero;
        }
    }
    //
    //     B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS)
    //
    if (rank < n) {
        Rormrz("Left", "Transpose", n, nrhs, rank, n - rank, a, lda, &work[(mn + 1) - 1], b, ldb, &work[(2 * mn + 1) - 1], lwork - 2 * mn, info);
    }
    //
    //     workspace: 2*MN+NRHS.
    //
    //     B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
    //
    for (j = 1; j <= nrhs; j = j + 1) {
        for (i = 1; i <= n; i = i + 1) {
            work[jpvt[i - 1] - 1] = b[(i - 1) + (j - 1) * ldb];
        }
        Rcopy(n, &work[1 - 1], 1, &b[(j - 1) * ldb], 1);
    }
    //
    //     workspace: N.
    //
    //     Undo scaling
    //
    if (iascl == 1) {
        Rlascl("G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
        Rlascl("U", 0, 0, smlnum, anrm, rank, rank, a, lda, info);
    } else if (iascl == 2) {
        Rlascl("G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
        Rlascl("U", 0, 0, bignum, anrm, rank, rank, a, lda, info);
    }
    if (ibscl == 1) {
        Rlascl("G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
    } else if (ibscl == 2) {
        Rlascl("G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
    }
//
statement_70:
    work[1 - 1] = lwkopt;
    //
    //     End of Rgelsy
    //
}
