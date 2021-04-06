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

void Rgeqrt(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *a, INTEGER const lda, REAL *t, INTEGER const ldt, REAL *work, INTEGER &info) {
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
    } else if (nb < 1 || (nb > min(m, n) && min(m, n) > 0)) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldt < nb) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Rgeqrt", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    INTEGER k = min(m, n);
    if (k == 0) {
        return;
    }
    //
    //     Blocked loop of length K
    //
    INTEGER i = 0;
    INTEGER ib = 0;
    const bool use_recursive_qr = true;
    INTEGER iinfo = 0;
    for (i = 1; i <= k; i = i + nb) {
        ib = min(k - i + 1, nb);
        //
        //     Compute the QR factorization of the current block A(I:M,I:I+IB-1)
        //
        if (use_recursive_qr) {
            Rgeqrt3(m - i + 1, ib, &a[(i - 1) + (i - 1) * lda], lda, &t[(i - 1) * ldt], ldt, iinfo);
        } else {
            Rgeqrt2(m - i + 1, ib, &a[(i - 1) + (i - 1) * lda], lda, &t[(i - 1) * ldt], ldt, iinfo);
        }
        if (i + ib <= n) {
            //
            //     Update by applying H**T to A(I:M,I+IB:N) from the left
            //
            Rlarfb("L", "T", "F", "C", m - i + 1, n - i - ib + 1, ib, &a[(i - 1) + (i - 1) * lda], lda, &t[(i - 1) * ldt], ldt, &a[(i - 1) + ((i + ib) - 1) * lda], lda, work, n - i - ib + 1);
        }
    }
    //
    //     End of Rgeqrt
    //
}
