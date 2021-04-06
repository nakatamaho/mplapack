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

void Rgelq(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *t, INTEGER const tsize, REAL *work, INTEGER const lwork, INTEGER &info) {
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
        mb = iMlaenv(1, "Rgelq ", " ", m, n, 1, -1);
        nb = iMlaenv(1, "Rgelq ", " ", m, n, 2, -1);
    } else {
        mb = 1;
        nb = n;
    }
    if (mb > min(m, n) || mb < 1) {
        mb = 1;
    }
    if (nb > n || nb <= m) {
        nb = n;
    }
    INTEGER mINTEGERsz = m + 5;
    INTEGER nblcks = 0;
    if (nb > m && n > m) {
        if (mod(n - m, nb - m) == 0) {
            nblcks = (n - m) / (nb - m);
        } else {
            nblcks = (n - m) / (nb - m) + 1;
        }
    } else {
        nblcks = 1;
    }
    //
    //     Determine if the workspace size satisfies minimal size
    //
    INTEGER lwmin = 0;
    INTEGER lwopt = 0;
    if ((n <= m) || (nb <= m) || (nb >= n)) {
        lwmin = max((INTEGER)1, n);
        lwopt = max((INTEGER)1, mb * n);
    } else {
        lwmin = max((INTEGER)1, m);
        lwopt = max((INTEGER)1, mb * m);
    }
    bool lminws = false;
    if ((tsize < max((INTEGER)1, mb * m * nblcks + 5) || lwork < lwopt) && (lwork >= lwmin) && (tsize >= mINTEGERsz) && (!lquery)) {
        if (tsize < max((INTEGER)1, mb * m * nblcks + 5)) {
            lminws = true;
            mb = 1;
            nb = n;
        }
        if (lwork < lwopt) {
            lminws = true;
            mb = 1;
        }
    }
    INTEGER lwreq = 0;
    if ((n <= m) || (nb <= m) || (nb >= n)) {
        lwreq = max((INTEGER)1, mb * n);
    } else {
        lwreq = max((INTEGER)1, mb * m);
    }
    //
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    } else if (tsize < max((INTEGER)1, mb * m * nblcks + 5) && (!lquery) && (!lminws)) {
        info = -6;
    } else if ((lwork < lwreq) && (!lquery) && (!lminws)) {
        info = -8;
    }
    //
    if (info == 0) {
        if (mINTEGER) {
            t[1 - 1] = mINTEGERsz;
        } else {
            t[1 - 1] = mb * m * nblcks + 5;
        }
        t[2 - 1] = mb;
        t[3 - 1] = nb;
        if (minw) {
            work[1 - 1] = lwmin;
        } else {
            work[1 - 1] = lwreq;
        }
    }
    if (info != 0) {
        Mxerbla("Rgelq", -info);
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
    //     The LQ Decomposition
    //
    if ((n <= m) || (nb <= m) || (nb >= n)) {
        Rgelqt(m, n, mb, a, lda, &t[6 - 1], mb, work, info);
    } else {
        Rlaswlq(m, n, mb, nb, a, lda, &t[6 - 1], mb, work, lwork, info);
    }
    //
    work[1 - 1] = lwreq;
    //
    //     End of Rgelq
    //
}
