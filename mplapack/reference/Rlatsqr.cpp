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

void Rlatsqr(INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *a, INTEGER const lda, REAL *t, INTEGER const ldt, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. EXTERNAL FUNCTIONS ..
    //     .. EXTERNAL SUBROUTINES ..
    //     .. INTRINSIC FUNCTIONS ..
    //     ..
    //     .. EXECUTABLE STATEMENTS ..
    //
    //     TEST THE INPUT ARGUMENTS
    //
    info = 0;
    //
    bool lquery = (lwork == -1);
    //
    if (m < 0) {
        info = -1;
    } else if (n < 0 || m < n) {
        info = -2;
    } else if (mb <= n) {
        info = -3;
    } else if (nb < 1 || (nb > n && n > 0)) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldt < nb) {
        info = -8;
    } else if (lwork < (n * nb) && (!lquery)) {
        info = -10;
    }
    if (info == 0) {
        work[1 - 1] = nb * n;
    }
    if (info != 0) {
        Mxerbla("Rlatsqr", -info);
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
    if ((mb <= n) || (mb >= m)) {
        Rgeqrt(m, n, nb, a, lda, t, ldt, work, info);
        return;
    }
    //
    INTEGER kk = mod((m - n), (mb - n));
    INTEGER ii = m - kk + 1;
    //
    //      Compute the QR factorization of the first block A(1:MB,1:N)
    //
    Rgeqrt(mb, n, nb, &a[(1 - 1)], lda, t, ldt, work, info);
    //
    INTEGER ctr = 1;
    INTEGER i = 0;
    for (i = mb + 1; i <= ii - mb + n; i = i + (mb - n)) {
        //
        //      Compute the QR factorization of the current block A(I:I+MB-N,1:N)
        //
        Rtpqrt(mb - n, n, 0, nb, &a[(1 - 1)], lda, &a[(i - 1)], lda, &t[((ctr * n + 1) - 1) * ldt], ldt, work, info);
        ctr++;
    }
    //
    //      Compute the QR factorization of the last block A(II:M,1:N)
    //
    if (ii <= m) {
        Rtpqrt(kk, n, 0, nb, &a[(1 - 1)], lda, &a[(ii - 1)], lda, &t[((ctr * n + 1) - 1) * ldt], ldt, work, info);
    }
    //
    work[1 - 1] = n * nb;
    //
    //     End of Rlatsqr
    //
}
