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

void Rormtr(const char *side, const char *uplo, const char *trans, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *tau, REAL *c, INTEGER const ldc, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    bool left = Mlsame(side, "L");
    bool upper = Mlsame(uplo, "U");
    bool lquery = (lwork == -1);
    //
    //
    INTEGER nq = 0;
    INTEGER nw = 0;
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if (!left && !Mlsame(side, "R")) {
        info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "T")) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, nq)) {
        info = -7;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -10;
    } else if (lwork < max((INTEGER)1, nw) && !lquery) {
        info = -12;
    }
    //
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        if (upper) {
            if (left) {
                nb = iMlaenv(1, "Rormql", side + trans, m - 1, n, m - 1, -1);
            } else {
                nb = iMlaenv(1, "Rormql", side + trans, m, n - 1, n - 1, -1);
            }
        } else {
            if (left) {
                nb = iMlaenv(1, "Rormqr", side + trans, m - 1, n, m - 1, -1);
            } else {
                nb = iMlaenv(1, "Rormqr", side + trans, m, n - 1, n - 1, -1);
            }
        }
        lwkopt = max((INTEGER)1, nw) * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rormtr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0 || nq == 1) {
        work[1 - 1] = 1;
        return;
    }
    //
    INTEGER mi = 0;
    INTEGER ni = 0;
    if (left) {
        mi = m - 1;
        ni = n;
    } else {
        mi = m;
        ni = n - 1;
    }
    //
    INTEGER iinfo = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    if (upper) {
        //
        //        Q was determined by a call to Rsytrd with UPLO = 'U'
        //
        Rormql(side, trans, mi, ni, nq - 1, &a[(2 - 1) * lda], lda, tau, c, ldc, work, lwork, iinfo);
    } else {
        //
        //        Q was determined by a call to Rsytrd with UPLO = 'L'
        //
        if (left) {
            i1 = 2;
            i2 = 1;
        } else {
            i1 = 1;
            i2 = 2;
        }
        Rormqr(side, trans, mi, ni, nq - 1, &a[(2 - 1)], lda, tau, &c[(i1 - 1) + (i2 - 1) * ldc], ldc, work, lwork, iinfo);
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Rormtr
    //
}
