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

void Rlaswlq(INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *a, INTEGER const lda, REAL *t, INTEGER const ldt, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    } else if (n < 0 || n < m) {
        info = -2;
    } else if (mb < 1 || (mb > m && m > 0)) {
        info = -3;
    } else if (nb <= m) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldt < mb) {
        info = -8;
    } else if ((lwork < m * mb) && (!lquery)) {
        info = -10;
    }
    if (info == 0) {
        work[1 - 1] = mb * m;
    }
    //
    if (info != 0) {
        Mxerbla("Rlaswlq", -info);
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
    if ((m >= n) || (nb <= m) || (nb >= n)) {
        Rgelqt(m, n, mb, a, lda, t, ldt, work, info);
        return;
    }
    //
    INTEGER kk = mod((n - m), (nb - m));
    INTEGER ii = n - kk + 1;
    //
    //      Compute the LQ factorization of the first block A(1:M,1:NB)
    //
    Rgelqt(m, nb, mb, &a[(1 - 1)], lda, t, ldt, work, info);
    INTEGER ctr = 1;
    //
    INTEGER i = 0;
    for (i = nb + 1; i <= ii - nb + m; i = i + (nb - m)) {
        //
        //      Compute the QR factorization of the current block A(1:M,I:I+NB-M)
        //
        Rtplqt(m, nb - m, 0, mb, &a[(1 - 1)], lda, &a[(i - 1) * lda], lda, &t[((ctr * m + 1) - 1) * ldt], ldt, work, info);
        ctr++;
    }
    //
    //     Compute the QR factorization of the last block A(1:M,II:N)
    //
    if (ii <= n) {
        Rtplqt(m, kk, 0, mb, &a[(1 - 1)], lda, &a[(ii - 1) * lda], lda, &t[((ctr * m + 1) - 1) * ldt], ldt, work, info);
    }
    //
    work[1 - 1] = m * mb;
    //
    //     End of Rlaswlq
    //
}
