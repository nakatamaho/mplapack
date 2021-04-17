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

void Rgeqp3(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *jpvt, REAL *tau, REAL *work, INTEGER const lwork, INTEGER &info) {
    bool lquery = false;
    INTEGER minmn = 0;
    INTEGER iws = 0;
    INTEGER lwkopt = 0;
    const INTEGER inb = 1;
    INTEGER nb = 0;
    INTEGER nfxd = 0;
    INTEGER j = 0;
    INTEGER na = 0;
    INTEGER sm = 0;
    INTEGER sn = 0;
    INTEGER sminmn = 0;
    INTEGER nbmin = 0;
    INTEGER nx = 0;
    const INTEGER ixover = 3;
    INTEGER minws = 0;
    const INTEGER inbmin = 2;
    INTEGER topbmn = 0;
    INTEGER jb = 0;
    INTEGER fjb = 0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test input arguments
    //  ====================
    //
    info = 0;
    lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    //
    if (info == 0) {
        minmn = min(m, n);
        if (minmn == 0) {
            iws = 1;
            lwkopt = 1;
        } else {
            iws = 3 * n + 1;
            nb = iMlaenv(inb, "Rgeqrf", " ", m, n, -1, -1);
            lwkopt = 2 * n + (n + 1) * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if ((lwork < iws) && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgeqp3", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Move initial columns up front.
    //
    nfxd = 1;
    for (j = 1; j <= n; j = j + 1) {
        if (jpvt[j - 1] != 0) {
            if (j != nfxd) {
                Rswap(m, &a[(j - 1) * lda], 1, &a[(nfxd - 1) * lda], 1);
                jpvt[j - 1] = jpvt[nfxd - 1];
                jpvt[nfxd - 1] = j;
            } else {
                jpvt[j - 1] = j;
            }
            nfxd++;
        } else {
            jpvt[j - 1] = j;
        }
    }
    nfxd = nfxd - 1;
    //
    //     Factorize fixed columns
    //  =======================
    //
    //     Compute the QR factorization of fixed columns and update
    //     remaining columns.
    //
    if (nfxd > 0) {
        na = min(m, nfxd);
        // CC      CALL Rgeqr2( M, NA, A, LDA, TAU, WORK, INFO )
        Rgeqrf(m, na, a, lda, tau, work, lwork, info);
        iws = max(iws, castINTEGER(work[1 - 1]));
        if (na < n) {
            // CC         CALL Rorm2r( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
            // CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            Rormqr("Left", "Transpose", m, n - na, na, a, lda, tau, &a[((na + 1) - 1) * lda], lda, work, lwork, info);
            iws = max(iws, castINTEGER(work[1 - 1]));
        }
    }
    //
    //     Factorize free columns
    //  ======================
    //
    if (nfxd < minmn) {
        //
        sm = m - nfxd;
        sn = n - nfxd;
        sminmn = minmn - nfxd;
        //
        //        Determine the block size.
        //
        nb = iMlaenv(inb, "Rgeqrf", " ", sm, sn, -1, -1);
        nbmin = 2;
        nx = 0;
        //
        if ((nb > 1) && (nb < sminmn)) {
            //
            //           Determine when to cross over from blocked to unblocked code.
            //
            nx = max((INTEGER)0, iMlaenv(ixover, "Rgeqrf", " ", sm, sn, -1, -1));
            //
            if (nx < sminmn) {
                //
                //              Determine if workspace is large enough for blocked code.
                //
                minws = 2 * sn + (sn + 1) * nb;
                iws = max(iws, minws);
                if (lwork < minws) {
                    //
                    //                 Not enough workspace to use optimal NB: Reduce NB and
                    //                 determine the minimum value of NB.
                    //
                    nb = (lwork - 2 * sn) / (sn + 1);
                    nbmin = max((INTEGER)2, iMlaenv(inbmin, "Rgeqrf", " ", sm, sn, -1, -1));
                    //
                }
            }
        }
        //
        //        Initialize partial column norms. The first N elements of work
        //        store the exact column norms.
        //
        for (j = nfxd + 1; j <= n; j = j + 1) {
            work[j - 1] = Rnrm2(sm, &a[((nfxd + 1) - 1) + (j - 1) * lda], 1);
            work[(n + j) - 1] = work[j - 1];
        }
        //
        if ((nb >= nbmin) && (nb < sminmn) && (nx < sminmn)) {
            //
            //           Use blocked code initially.
            //
            j = nfxd + 1;
            //
            //           Compute factorization: while loop.
            //
            topbmn = minmn - nx;
        statement_30:
            if (j <= topbmn) {
                jb = min(nb, topbmn - j + 1);
                //
                //              Factorize JB columns among columns J:N.
                //
                Rlaqps(m, n - j + 1, j - 1, jb, fjb, &a[(j - 1) * lda], lda, &jpvt[j - 1], &tau[j - 1], &work[j - 1], &work[(n + j) - 1], &work[(2 * n + 1) - 1], &work[(2 * n + jb + 1) - 1], n - j + 1);
                //
                j += fjb;
                goto statement_30;
            }
        } else {
            j = nfxd + 1;
        }
        //
        //        Use unblocked code to factor the last or only block.
        //
        if (j <= minmn) {
            Rlaqp2(m, n - j + 1, j - 1, &a[(j - 1) * lda], lda, &jpvt[j - 1], &tau[j - 1], &work[j - 1], &work[(n + j) - 1], &work[(2 * n + 1) - 1]);
        }
        //
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Rgeqp3
    //
}
