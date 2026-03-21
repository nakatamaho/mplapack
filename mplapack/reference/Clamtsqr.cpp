/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZLAMTSQR.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Clamtsqr(const char *side, const char *trans, INTEGER const m, INTEGER const n, INTEGER const k, INTEGER const mb, INTEGER const nb, COMPLEX *a, INTEGER const lda, COMPLEX *t, INTEGER const ldt, COMPLEX *c, INTEGER const ldc, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
    //
    // Test the input arguments
    //
    info = 0;
    bool lquery = (lwork == -1);
    bool notran = Mlsame(trans, "N");
    bool tran = Mlsame(trans, "C");
    bool left = Mlsame(side, "L");
    bool right = Mlsame(side, "R");
    INTEGER lw = 0;
    INTEGER q = 0;
    if (left) {
        lw = n * nb;
        q = m;
    } else {
        lw = m * nb;
        q = n;
    }
    //
    INTEGER minmnk = min(m, n, k);
    INTEGER lwmin = 0;
    if (minmnk == 0) {
        lwmin = 1;
    } else {
        lwmin = max((INTEGER)1, lw);
    }
    //
    if (!left && !right) {
        info = -1;
    } else if (!tran && !notran) {
        info = -2;
    } else if (m < k) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0) {
        info = -5;
    } else if (k < nb || nb < 1) {
        info = -7;
    } else if (lda < max((INTEGER)1, q)) {
        info = -9;
    } else if (ldt < max((INTEGER)1, nb)) {
        info = -11;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -13;
    } else if (lwork < lwmin && (!lquery)) {
        info = -15;
    }
    //
    if (info == 0) {
        work[1 - 1] = lwmin;
    }
    //
    if (info != 0) {
        Mxerbla("Clamtsqr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    // Quick return if possible
    //
    if (minmnk == 0) {
        return;
    }
    //
    // Determine the block size if it is tall skinny or short and wide
    //
    if ((mb <= k) || (mb >= max(m, n, k))) {
        Cgemqrt(side, trans, m, n, k, nb, a, lda, t, ldt, c, ldc, work, info);
        return;
    }
    //
    INTEGER kk = 0;
    INTEGER ctr = 0;
    INTEGER ii = 0;
    INTEGER i = 0;
    if (left && notran) {
        //
        // Multiply Q to the last block of C
        //
        kk = mod((m - k), (mb - k));
        ctr = (m - k) / (mb - k);
        if (kk > 0) {
            ii = m - kk + 1;
            Ctpmqrt("L", "N", kk, n, k, 0, nb, &a[(ii - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(ii - 1)], ldc, work, info);
        } else {
            ii = m + 1;
        }
        //
        for (i = ii - (mb - k); i >= mb + 1; i = i - (mb - k)) {
            //
            // Multiply Q to the current block of C (I:I+MB,1:N)
            //
            ctr = ctr - 1;
            Ctpmqrt("L", "N", mb - k, n, k, 0, nb, &a[(i - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(i - 1)], ldc, work, info);
            //
        }
        //
        // Multiply Q to the first block of C (1:MB,1:N)
        //
        Cgemqrt("L", "N", mb, n, k, nb, &a[0], lda, t, ldt, &c[0], ldc, work, info);
        //
    } else if (left && tran) {
        //
        // Multiply Q to the first block of C
        //
        kk = mod((m - k), (mb - k));
        ii = m - kk + 1;
        ctr = 1;
        Cgemqrt("L", "C", mb, n, k, nb, &a[0], lda, t, ldt, &c[0], ldc, work, info);
        //
        for (i = mb + 1; i <= ii - mb + k; i = i + (mb - k)) {
            //
            // Multiply Q to the current block of C (I:I+MB,1:N)
            //
            Ctpmqrt("L", "C", mb - k, n, k, 0, nb, &a[(i - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(i - 1)], ldc, work, info);
            ctr++;
            //
        }
        if (ii <= m) {
            //
            // Multiply Q to the last block of C
            //
            Ctpmqrt("L", "C", kk, n, k, 0, nb, &a[(ii - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(ii - 1)], ldc, work, info);
            //
        }
        //
    } else if (right && tran) {
        //
        // Multiply Q to the last block of C
        //
        kk = mod((n - k), (mb - k));
        ctr = (n - k) / (mb - k);
        if (kk > 0) {
            ii = n - kk + 1;
            Ctpmqrt("R", "C", m, kk, k, 0, nb, &a[(ii - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(ii - 1) * ldc], ldc, work, info);
        } else {
            ii = n + 1;
        }
        //
        for (i = ii - (mb - k); i >= mb + 1; i = i - (mb - k)) {
            //
            // Multiply Q to the current block of C (1:M,I:I+MB)
            //
            ctr = ctr - 1;
            Ctpmqrt("R", "C", m, mb - k, k, 0, nb, &a[(i - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(i - 1) * ldc], ldc, work, info);
            //
        }
        //
        // Multiply Q to the first block of C (1:M,1:MB)
        //
        Cgemqrt("R", "C", m, mb, k, nb, &a[0], lda, t, ldt, &c[0], ldc, work, info);
        //
    } else if (right && notran) {
        //
        // Multiply Q to the first block of C
        //
        kk = mod((n - k), (mb - k));
        ii = n - kk + 1;
        ctr = 1;
        Cgemqrt("R", "N", m, mb, k, nb, &a[0], lda, t, ldt, &c[0], ldc, work, info);
        //
        for (i = mb + 1; i <= ii - mb + k; i = i + (mb - k)) {
            //
            // Multiply Q to the current block of C (1:M,I:I+MB)
            //
            Ctpmqrt("R", "N", m, mb - k, k, 0, nb, &a[(i - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(i - 1) * ldc], ldc, work, info);
            ctr++;
            //
        }
        if (ii <= n) {
            //
            // Multiply Q to the last block of C
            //
            Ctpmqrt("R", "N", m, kk, k, 0, nb, &a[(ii - 1)], lda, &t[((ctr * k + 1) - 1) * ldt], ldt, &c[0], ldc, &c[(ii - 1) * ldc], ldc, work, info);
            //
        }
        //
    }
    //
    work[1 - 1] = lwmin;
    //
    // End of Clamtsqr
    //
}
