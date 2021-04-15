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

void Rgetsls(const char *trans, INTEGER const m, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *work, INTEGER const lwork, INTEGER &info) {
    INTEGER minmn = 0;
    INTEGER maxmn = 0;
    INTEGER mnk = 0;
    bool tran = false;
    bool lquery = false;
    REAL tq[5];
    REAL workq[1];
    INTEGER info2 = 0;
    INTEGER tszo = 0;
    INTEGER lwo = 0;
    INTEGER tszm = 0;
    INTEGER lwm = 0;
    INTEGER wsizeo = 0;
    INTEGER wsizem = 0;
    INTEGER lw1 = 0;
    INTEGER lw2 = 0;
    const REAL zero = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    INTEGER iascl = 0;
    INTEGER brow = 0;
    REAL bnrm = 0.0;
    INTEGER ibscl = 0;
    INTEGER scllen = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    //
    //  -- LAPACK driver routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //
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
    //     Test the input arguments.
    //
    info = 0;
    minmn = min(m, n);
    maxmn = max(m, n);
    mnk = max(minmn, nrhs);
    tran = Mlsame(trans, "T");
    //
    lquery = (lwork == -1 || lwork == -2);
    if (!(Mlsame(trans, "N") || Mlsame(trans, "T"))) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (nrhs < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (ldb < max({(INTEGER)1, m, n})) {
        info = -8;
    }
    //
    if (info == 0) {
        //
        //     Determine the block size and minimum LWORK
        //
        if (m >= n) {
            Rgeqr(m, n, a, lda, tq, -1, workq, -1, info2);
            tszo = castINTEGER(tq[1 - 1]);
            lwo = castINTEGER(workq[1 - 1]);
            Rgemqr("L", trans, m, nrhs, n, a, lda, tq, tszo, b, ldb, workq, -1, info2);
            lwo = max(lwo, castINTEGER(workq[1 - 1]));
            Rgeqr(m, n, a, lda, tq, -2, workq, -2, info2);
            tszm = castINTEGER(tq[1 - 1]);
            lwm = castINTEGER(workq[1 - 1]);
            Rgemqr("L", trans, m, nrhs, n, a, lda, tq, tszm, b, ldb, workq, -1, info2);
            lwm = max(lwm, castINTEGER(workq[1 - 1]));
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        } else {
            Rgelq(m, n, a, lda, tq, -1, workq, -1, info2);
            tszo = castINTEGER(tq[1 - 1]);
            lwo = castINTEGER(workq[1 - 1]);
            Rgemlq("L", trans, n, nrhs, m, a, lda, tq, tszo, b, ldb, workq, -1, info2);
            lwo = max(lwo, castINTEGER(workq[1 - 1]));
            Rgelq(m, n, a, lda, tq, -2, workq, -2, info2);
            tszm = castINTEGER(tq[1 - 1]);
            lwm = castINTEGER(workq[1 - 1]);
            Rgemlq("L", trans, n, nrhs, m, a, lda, tq, tszm, b, ldb, workq, -1, info2);
            lwm = max(lwm, castINTEGER(workq[1 - 1]));
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        }
        //
        if ((lwork < wsizem) && (!lquery)) {
            info = -10;
        }
        //
    }
    //
    if (info != 0) {
        Mxerbla("Rgetsls", -info);
        work[1 - 1] = castREAL(wsizeo);
        return;
    }
    if (lquery) {
        if (lwork == -1) {
            work[1 - 1] = castREAL(wsizeo);
        }
        if (lwork == -2) {
            work[1 - 1] = castREAL(wsizem);
        }
        return;
    }
    if (lwork < wsizeo) {
        lw1 = tszm;
        lw2 = lwm;
    } else {
        lw1 = tszo;
        lw2 = lwo;
    }
    //
    //     Quick return if possible
    //
    if (min({m, n, nrhs}) == 0) {
        Rlaset("FULL", max(m, n), nrhs, zero, zero, b, ldb);
        return;
    }
    //
    //     Get machine parameters
    //
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Scale A, B if max element outside range [SMLNUM,BIGNUM]
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
        Rlaset("F", maxmn, nrhs, zero, zero, b, ldb);
        goto statement_50;
    }
    //
    brow = m;
    if (tran) {
        brow = n;
    }
    bnrm = Rlange("M", brow, nrhs, b, ldb, work);
    ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Rlascl("G", 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Rlascl("G", 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info);
        ibscl = 2;
    }
    //
    if (m >= n) {
        //
        //        compute QR factorization of A
        //
        Rgeqr(m, n, a, lda, &work[(lw2 + 1) - 1], lw1, &work[1 - 1], lw2, info);
        if (!tran) {
            //
            //           Least-Squares Problem min || A * X - B ||
            //
            //           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
            //
            Rgemqr("L", "T", m, nrhs, n, a, lda, &work[(lw2 + 1) - 1], lw1, b, ldb, &work[1 - 1], lw2, info);
            //
            //           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
            //
            Rtrtrs("U", "N", "N", n, nrhs, a, lda, b, ldb, info);
            if (info > 0) {
                return;
            }
            scllen = n;
        } else {
            //
            //           Overdetermined system of equations A**T * X = B
            //
            //           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
            //
            Rtrtrs("U", "T", "N", n, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            //           B(N+1:M,1:NRHS) = ZERO
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                for (i = n + 1; i <= m; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = zero;
                }
            }
            //
            //           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
            //
            Rgemqr("L", "N", m, nrhs, n, a, lda, &work[(lw2 + 1) - 1], lw1, b, ldb, &work[1 - 1], lw2, info);
            //
            scllen = m;
            //
        }
        //
    } else {
        //
        //        Compute LQ factorization of A
        //
        Rgelq(m, n, a, lda, &work[(lw2 + 1) - 1], lw1, &work[1 - 1], lw2, info);
        //
        //        workspace at least M, optimally M*NB.
        //
        if (!tran) {
            //
            //           underdetermined system of equations A * X = B
            //
            //           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
            //
            Rtrtrs("L", "N", "N", m, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            //           B(M+1:N,1:NRHS) = 0
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                for (i = m + 1; i <= n; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = zero;
                }
            }
            //
            //           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
            //
            Rgemlq("L", "T", n, nrhs, m, a, lda, &work[(lw2 + 1) - 1], lw1, b, ldb, &work[1 - 1], lw2, info);
            //
            //           workspace at least NRHS, optimally NRHS*NB
            //
            scllen = n;
            //
        } else {
            //
            //           overdetermined system min || A**T * X - B ||
            //
            //           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
            //
            Rgemlq("L", "N", n, nrhs, m, a, lda, &work[(lw2 + 1) - 1], lw1, b, ldb, &work[1 - 1], lw2, info);
            //
            //           workspace at least NRHS, optimally NRHS*NB
            //
            //           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
            //
            Rtrtrs("Lower", "Transpose", "Non-unit", m, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            scllen = m;
            //
        }
        //
    }
    //
    //     Undo scaling
    //
    if (iascl == 1) {
        Rlascl("G", 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info);
    } else if (iascl == 2) {
        Rlascl("G", 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info);
    }
    if (ibscl == 1) {
        Rlascl("G", 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info);
    } else if (ibscl == 2) {
        Rlascl("G", 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info);
    }
//
statement_50:
    work[1 - 1] = castREAL(tszo + lwo);
    //
    //     End of Rgetsls
    //
}
