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

void Rlamtsqr(const char *side, const char *trans, INTEGER const &m, INTEGER const &n, INTEGER const &k, INTEGER const &mb, INTEGER const &nb, REAL *a, INTEGER const &lda, REAL *t, INTEGER const &ldt, REAL *c, INTEGER const &ldc, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
    //     .. External Functions ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    bool lquery = lwork < 0;
    bool notran = Mlsame(trans, "N");
    bool tran = Mlsame(trans, "T");
    bool left = Mlsame(side, "L");
    bool right = Mlsame(side, "R");
    INTEGER lw = 0;
    if (left) {
        lw = n * nb;
    } else {
        lw = mb * nb;
    }
    //
    info = 0;
    if (!left && !right) {
        info = -1;
    } else if (!tran && !notran) {
        info = -2;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, k)) {
        info = -9;
    } else if (ldt < max((INTEGER)1, nb)) {
        info = -11;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -13;
    } else if ((lwork < max((INTEGER)1, lw)) && (!lquery)) {
        info = -15;
    }
    //
    //     Determine the block size if it is tall skinny or short and wide
    //
    if (info == 0) {
        work[1 - 1] = lw;
    }
    //
    if (info != 0) {
        Mxerbla("Rlamtsqr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n, k) == 0) {
        return;
    }
    //
    if ((mb <= k) || (mb >= max(m, n, k))) {
        Rgemqrt(side, trans, m, n, k, nb, a, lda, t, ldt, c, ldc, work, info);
        return;
    }
    //
    INTEGER kk = 0;
    INTEGER ctr = 0;
    INTEGER ii = 0;
    INTEGER i = 0;
    if (left && notran) {
        //
        //         Multiply Q to the last block of C
        //
        kk = mod((m - k), (mb - k));
        ctr = (m - k) / (mb - k);
        if (kk > 0) {
            ii = m - kk + 1;
            Rtpmqrt("L", "N", kk, n, k, 0, nb, a[(ii - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(ii - 1)], ldc, work, info);
        } else {
            ii = m + 1;
        }
        //
        for (i = ii - (mb - k); i <= mb + 1; i = i + -(mb - k)) {
            //
            //         Multiply Q to the current block of C (I:I+MB,1:N)
            //
            ctr = ctr - 1;
            Rtpmqrt("L", "N", mb - k, n, k, 0, nb, a[(i - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(i - 1)], ldc, work, info);
            //
        }
        //
        //         Multiply Q to the first block of C (1:MB,1:N)
        //
        Rgemqrt("L", "N", mb, n, k, nb, a[(1 - 1)], lda, t, ldt, c[(1 - 1)], ldc, work, info);
        //
    } else if (left && tran) {
        //
        //         Multiply Q to the first block of C
        //
        kk = mod((m - k), (mb - k));
        ii = m - kk + 1;
        ctr = 1;
        Rgemqrt("L", "T", mb, n, k, nb, a[(1 - 1)], lda, t, ldt, c[(1 - 1)], ldc, work, info);
        //
        for (i = mb + 1; i <= ii - mb + k; i = i + (mb - k)) {
            //
            //         Multiply Q to the current block of C (I:I+MB,1:N)
            //
            Rtpmqrt("L", "T", mb - k, n, k, 0, nb, a[(i - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(i - 1)], ldc, work, info);
            ctr++;
            //
        }
        if (ii <= m) {
            //
            //         Multiply Q to the last block of C
            //
            Rtpmqrt("L", "T", kk, n, k, 0, nb, a[(ii - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(ii - 1)], ldc, work, info);
            //
        }
        //
    } else if (right && tran) {
        //
        //         Multiply Q to the last block of C
        //
        kk = mod((n - k), (mb - k));
        ctr = (n - k) / (mb - k);
        if (kk > 0) {
            ii = n - kk + 1;
            Rtpmqrt("R", "T", m, kk, k, 0, nb, a[(ii - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(ii - 1) * ldc], ldc, work, info);
        } else {
            ii = n + 1;
        }
        //
        for (i = ii - (mb - k); i <= mb + 1; i = i + -(mb - k)) {
            //
            //         Multiply Q to the current block of C (1:M,I:I+MB)
            //
            ctr = ctr - 1;
            Rtpmqrt("R", "T", m, mb - k, k, 0, nb, a[(i - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(i - 1) * ldc], ldc, work, info);
            //
        }
        //
        //         Multiply Q to the first block of C (1:M,1:MB)
        //
        Rgemqrt("R", "T", m, mb, k, nb, a[(1 - 1)], lda, t, ldt, c[(1 - 1)], ldc, work, info);
        //
    } else if (right && notran) {
        //
        //         Multiply Q to the first block of C
        //
        kk = mod((n - k), (mb - k));
        ii = n - kk + 1;
        ctr = 1;
        Rgemqrt("R", "N", m, mb, k, nb, a[(1 - 1)], lda, t, ldt, c[(1 - 1)], ldc, work, info);
        //
        for (i = mb + 1; i <= ii - mb + k; i = i + (mb - k)) {
            //
            //         Multiply Q to the current block of C (1:M,I:I+MB)
            //
            Rtpmqrt("R", "N", m, mb - k, k, 0, nb, a[(i - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(i - 1) * ldc], ldc, work, info);
            ctr++;
            //
        }
        if (ii <= n) {
            //
            //         Multiply Q to the last block of C
            //
            Rtpmqrt("R", "N", m, kk, k, 0, nb, a[(ii - 1)], lda, t[((ctr * k + 1) - 1) * ldt], ldt, c[(1 - 1)], ldc, c[(ii - 1) * ldc], ldc, work, info);
            //
        }
        //
    }
    //
    work[1 - 1] = lw;
    //
    //     End of Rlamtsqr
    //
}
