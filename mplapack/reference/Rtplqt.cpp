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

void Rtplqt(INTEGER const &m, INTEGER const &n, INTEGER const &l, INTEGER const &mb, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *t, INTEGER const &ldt, REAL *work, INTEGER &info) {
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
    // =====================================================================
    //
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (l < 0 || (l > min(m, n) && min(m, n) >= 0)) {
        info = -3;
    } else if (mb < 1 || (mb > m && m > 0)) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (ldb < max((INTEGER)1, m)) {
        info = -8;
    } else if (ldt < mb) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Rtplqt", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    INTEGER i = 0;
    INTEGER ib = 0;
    INTEGER nb = 0;
    INTEGER lb = 0;
    INTEGER iinfo = 0;
    for (i = 1; i <= m; i = i + mb) {
        //
        //     Compute the QR factorization of the current block
        //
        ib = min(m - i + 1, mb);
        nb = min(n - l + i + ib - 1, n);
        if (i >= l) {
            lb = 0;
        } else {
            lb = nb - n + l - i + 1;
        }
        //
        Rtplqt2(ib, nb, lb, a[(i - 1) + (i - 1) * lda], lda, b[(i - 1)], ldb, t[(i - 1) * ldt], ldt, iinfo);
        //
        //     Update by applying H**T to B(I+IB:M,:) from the right
        //
        if (i + ib <= m) {
            Rtprfb("R", "N", "F", "R", m - i - ib + 1, nb, ib, lb, b[(i - 1)], ldb, t[(i - 1) * ldt], ldt, a[((i + ib) - 1) + (i - 1) * lda], lda, b[((i + ib) - 1)], ldb, work, m - i - ib + 1);
        }
    }
    //
    //     End of Rtplqt
    //
}
