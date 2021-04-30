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

void Cggqrf(INTEGER const n, INTEGER const m, INTEGER const p, COMPLEX *a, INTEGER const lda, COMPLEX *taua, COMPLEX *b, INTEGER const ldb, COMPLEX *taub, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //     Test the input parameters
    //
    info = 0;
    INTEGER nb1 = iMlaenv(1, "Cgeqrf", " ", n, m, -1, -1);
    INTEGER nb2 = iMlaenv(1, "Cgerqf", " ", n, p, -1, -1);
    INTEGER nb3 = iMlaenv(1, "Cunmqr", " ", n, m, p, -1);
    INTEGER nb = max({nb1, nb2, nb3});
    INTEGER lwkopt = max({n, m, p}) * nb;
    work[1 - 1] = lwkopt;
    bool lquery = (lwork == -1);
    if (n < 0) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (p < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    } else if (lwork < max({(INTEGER)1, n, m, p}) && !lquery) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Cggqrf", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     QR factorization of N-by-M matrix A: A = Q*R
    //
    Cgeqrf(n, m, a, lda, taua, work, lwork, info);
    INTEGER lopt = castINTEGER(work[1 - 1].real());
    //
    //     Update B := Q**H*B.
    //
    Cunmqr("Left", "Conjugate Transpose", n, p, min(n, m), a, lda, taua, b, ldb, work, lwork, info);
    lopt = max(lopt, castINTEGER(work[1 - 1].real()));
    //
    //     RQ factorization of N-by-P matrix B: B = T*Z.
    //
    Cgerqf(n, p, b, ldb, taub, work, lwork, info);
    work[1 - 1] = max(lopt, castINTEGER(work[1 - 1].real()));
    //
    //     End of Cggqrf
    //
}
