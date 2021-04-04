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

void Rgeqr(INTEGER const &m, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *t, INTEGER const &tsize, REAL *work, INTEGER const &lwork, INTEGER &info) {
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    //
    bool lquery = (tsize == -1 || tsize == -2 || lwork == -1 || lwork == -2);
    //
    bool mINTEGER = false;
    bool minw = false;
    if (tsize == -2 || lwork == -2) {
        if (tsize != -1) {
            mINTEGER = true;
        }
        if (lwork != -1) {
            minw = true;
        }
    }
    //
    //     Determine the block size
    //
    INTEGER mb = 0;
    INTEGER nb = 0;
    if (min(m, n) > 0) {
        mb = iMlaenv[(("Rgeqr ") - 1) * ldiMlaenv];
        nb = iMlaenv[(("Rgeqr ") - 1) * ldiMlaenv];
    } else {
        mb = m;
        nb = 1;
    }
    if (mb > m || mb <= n) {
        mb = m;
    }
    if (nb > min(m, n) || nb < 1) {
        nb = 1;
    }
    INTEGER mINTEGERsz = n + 5;
    INTEGER nblcks = 0;
    if (mb > n && m > n) {
        if (mod(m - n, mb - n) == 0) {
            nblcks = (m - n) / (mb - n);
        } else {
            nblcks = (m - n) / (mb - n) + 1;
        }
    } else {
        nblcks = 1;
    }
    //
    //     Determine if the workspace size satisfies minimal size
    //
    bool lminws = false;
    if ((tsize < max((INTEGER)1, nb * n * nblcks + 5) || lwork < nb * n) && (lwork >= n) && (tsize >= mINTEGERsz) && (!lquery)) {
        if (tsize < max((INTEGER)1, nb * n * nblcks + 5)) {
            lminws = true;
            nb = 1;
            mb = m;
        }
        if (lwork < nb * n) {
            lminws = true;
            nb = 1;
        }
    }
    //
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    } else if (tsize < max((INTEGER)1, nb * n * nblcks + 5) && (!lquery) && (!lminws)) {
        info = -6;
    } else if ((lwork < max((INTEGER)1, n * nb)) && (!lquery) && (!lminws)) {
        info = -8;
    }
    //
    if (info == 0) {
        if (mINTEGER) {
            t[1 - 1] = mINTEGERsz;
        } else {
            t[1 - 1] = nb * n * nblcks + 5;
        }
        t[2 - 1] = mb;
        t[3 - 1] = nb;
        if (minw) {
            work[1 - 1] = max((INTEGER)1, n);
        } else {
            work[1 - 1] = max((INTEGER)1, nb * n);
        }
    }
    if (info != 0) {
        Mxerbla("Rgeqr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n) == 0) {
        return;
    }
    //
    //     The QR Decomposition
    //
    if ((m <= n) || (mb <= n) || (mb >= m)) {
        Rgeqrt(m, n, nb, a, lda, t[6 - 1], nb, work, info);
    } else {
        Rlatsqr(m, n, mb, nb, a, lda, t[6 - 1], nb, work, lwork, info);
    }
    //
    work[1 - 1] = max((INTEGER)1, nb * n);
    //
    //     End of Rgeqr
    //
}
