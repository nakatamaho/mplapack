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

void Cgemlqt(const char *side, const char *trans, INTEGER const &m, INTEGER const &n, INTEGER const &k, INTEGER const &mb, COMPLEX *v, INTEGER const &ldv, COMPLEX *t, INTEGER const &ldt, COMPLEX *c, INTEGER const &ldc, COMPLEX *work, INTEGER &info) {
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
    //     ..
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
    //     .. Test the input arguments ..
    //
    info = 0;
    bool left = Mlsame(side, "L");
    bool right = Mlsame(side, "R");
    bool tran = Mlsame(trans, "C");
    bool notran = Mlsame(trans, "N");
    //
    INTEGER ldwork = 0;
    if (left) {
        ldwork = max((INTEGER)1, n);
    } else if (right) {
        ldwork = max((INTEGER)1, m);
    }
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
    } else if (mb < 1 || (mb > k && k > 0)) {
        info = -6;
    } else if (ldv < max((INTEGER)1, k)) {
        info = -8;
    } else if (ldt < mb) {
        info = -10;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -12;
    }
    //
    if (info != 0) {
        Mxerbla("Cgemlqt", -info);
        return;
    }
    //
    //     .. Quick return if possible ..
    //
    if (m == 0 || n == 0 || k == 0) {
        return;
    }
    //
    INTEGER i = 0;
    INTEGER ib = 0;
    INTEGER kf = 0;
    if (left && notran) {
        //
        for (i = 1; i <= k; i = i + mb) {
            ib = min(mb, k - i + 1);
            Clarfb("L", "C", "F", "R", m - i + 1, n, ib, v[(i - 1) + (i - 1) * ldv], ldv, t[(i - 1) * ldt], ldt, c[(i - 1)], ldc, work, ldwork);
        }
        //
    } else if (right && tran) {
        //
        for (i = 1; i <= k; i = i + mb) {
            ib = min(mb, k - i + 1);
            Clarfb("R", "N", "F", "R", m, n - i + 1, ib, v[(i - 1) + (i - 1) * ldv], ldv, t[(i - 1) * ldt], ldt, c[(i - 1) * ldc], ldc, work, ldwork);
        }
        //
    } else if (left && tran) {
        //
        kf = ((k - 1) / mb) * mb + 1;
        for (i = kf; i <= 1; i = i + -mb) {
            ib = min(mb, k - i + 1);
            Clarfb("L", "N", "F", "R", m - i + 1, n, ib, v[(i - 1) + (i - 1) * ldv], ldv, t[(i - 1) * ldt], ldt, c[(i - 1)], ldc, work, ldwork);
        }
        //
    } else if (right && notran) {
        //
        kf = ((k - 1) / mb) * mb + 1;
        for (i = kf; i <= 1; i = i + -mb) {
            ib = min(mb, k - i + 1);
            Clarfb("R", "C", "F", "R", m, n - i + 1, ib, v[(i - 1) + (i - 1) * ldv], ldv, t[(i - 1) * ldt], ldt, c[(i - 1) * ldc], ldc, work, ldwork);
        }
        //
    }
    //
    //     End of Cgemlqt
    //
}
