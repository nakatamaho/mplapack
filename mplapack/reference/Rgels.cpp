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

void Rgels(const char *trans, INTEGER const m, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *work, INTEGER const lwork, INTEGER &info) {
    INTEGER mn = 0;
    bool lquery = false;
    bool tpsd = false;
    INTEGER nb = 0;
    INTEGER wsize = 0;
    const REAL zero = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL rwork[1];
    REAL anrm = 0.0;
    INTEGER iascl = 0;
    INTEGER brow = 0;
    REAL bnrm = 0.0;
    INTEGER ibscl = 0;
    INTEGER scllen = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    //
    //     Test the input arguments.
    //
    info = 0;
    mn = min(m, n);
    lquery = (lwork == -1);
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
    } else if (lwork < max((INTEGER)1, mn + max(mn, nrhs)) && !lquery) {
        info = -10;
    }
    //
    //     Figure out optimal block size
    //
    if (info == 0 || info == -10) {
        //
        tpsd = true;
        if (Mlsame(trans, "N")) {
            tpsd = false;
        }
        //
        if (m >= n) {
            nb = iMlaenv(1, "Rgeqrf", " ", m, n, -1, -1);
            if (tpsd) {
                nb = max(nb, iMlaenv(1, "Rormqr", "LN", m, nrhs, n, -1));
            } else {
                nb = max(nb, iMlaenv(1, "Rormqr", "LT", m, nrhs, n, -1));
            }
        } else {
            nb = iMlaenv(1, "Rgelqf", " ", m, n, -1, -1);
            if (tpsd) {
                nb = max(nb, iMlaenv(1, "Rormlq", "LT", n, nrhs, m, -1));
            } else {
                nb = max(nb, iMlaenv(1, "Rormlq", "LN", n, nrhs, m, -1));
            }
        }
        //
        wsize = max((INTEGER)1, mn + max(mn, nrhs) * nb);
        work[1 - 1] = castREAL(wsize);
        //
    }
    //
    if (info != 0) {
        Mxerbla("Rgels", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (min({m, n, nrhs}) == 0) {
        Rlaset("Full", max(m, n), nrhs, zero, zero, b, ldb);
        return;
    }
    //
    //     Get machine parameters
    //
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = one / smlnum;
    //
    //     Scale A, B if max element outside range [SMLNUM,BIGNUM]
    //
    anrm = Rlange("M", m, n, a, lda, rwork);
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
        goto statement_50;
    }
    //
    brow = m;
    if (tpsd) {
        brow = n;
    }
    bnrm = Rlange("M", brow, nrhs, b, ldb, rwork);
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
        Rgeqrf(m, n, a, lda, &work[1 - 1], &work[(mn + 1) - 1], lwork - mn, info);
        //
        //        workspace at least N, optimally N*NB
        //
        if (!tpsd) {
            //
            //           Least-Squares Problem min || A * X - B ||
            //
            //           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
            //
            Rormqr("Left", "Transpose", m, nrhs, n, a, lda, &work[1 - 1], b, ldb, &work[(mn + 1) - 1], lwork - mn, info);
            //
            //           workspace at least NRHS, optimally NRHS*NB
            //
            //           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
            //
            Rtrtrs("Upper", "No transpose", "Non-unit", n, nrhs, a, lda, b, ldb, info);
            //
            if (info > 0) {
                return;
            }
            //
            scllen = n;
            //
        } else {
            //
            //           Underdetermined system of equations A**T * X = B
            //
            //           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
            //
            Rtrtrs("Upper", "Transpose", "Non-unit", n, nrhs, a, lda, b, ldb, info);
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
            Rormqr("Left", "No transpose", m, nrhs, n, a, lda, &work[1 - 1], b, ldb, &work[(mn + 1) - 1], lwork - mn, info);
            //
            //           workspace at least NRHS, optimally NRHS*NB
            //
            scllen = m;
            //
        }
        //
    } else {
        //
        //        Compute LQ factorization of A
        //
        Rgelqf(m, n, a, lda, &work[1 - 1], &work[(mn + 1) - 1], lwork - mn, info);
        //
        //        workspace at least M, optimally M*NB.
        //
        if (!tpsd) {
            //
            //           underdetermined system of equations A * X = B
            //
            //           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
            //
            Rtrtrs("Lower", "No transpose", "Non-unit", m, nrhs, a, lda, b, ldb, info);
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
            Rormlq("Left", "Transpose", n, nrhs, m, a, lda, &work[1 - 1], b, ldb, &work[(mn + 1) - 1], lwork - mn, info);
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
            Rormlq("Left", "No transpose", n, nrhs, m, a, lda, &work[1 - 1], b, ldb, &work[(mn + 1) - 1], lwork - mn, info);
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
    work[1 - 1] = castREAL(wsize);
    //
    //     End of Rgels
    //
}
