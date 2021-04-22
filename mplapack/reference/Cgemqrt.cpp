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

void Cgemqrt(const char *side, const char *trans, INTEGER const m, INTEGER const n, INTEGER const k, INTEGER const nb, COMPLEX *v, INTEGER const ldv, COMPLEX *t, INTEGER const ldt, COMPLEX *c, INTEGER const ldc, COMPLEX *work, INTEGER &info) {
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
    INTEGER q = 0;
    if (left) {
        ldwork = max((INTEGER)1, n);
        q = m;
    } else if (right) {
        ldwork = max((INTEGER)1, m);
        q = n;
    }
    if (!left && !right) {
        info = -1;
    } else if (!tran && !notran) {
        info = -2;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0 || k > q) {
        info = -5;
    } else if (nb < 1 || (nb > k && k > 0)) {
        info = -6;
    } else if (ldv < max((INTEGER)1, q)) {
        info = -8;
    } else if (ldt < nb) {
        info = -10;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -12;
    }
    //
    if (info != 0) {
        Mxerbla("Cgemqrt", -info);
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
    if (left && tran) {
        //
        for (i = 1; i <= k; i = i + nb) {
            ib = min(nb, k - i + 1);
            Clarfb("L", "C", "F", "C", m - i + 1, n, ib, &v[(i - 1) + (i - 1) * ldv], ldv, &t[(i - 1) * ldt], ldt, &c[(i - 1)], ldc, work, ldwork);
        }
        //
    } else if (right && notran) {
        //
        for (i = 1; i <= k; i = i + nb) {
            ib = min(nb, k - i + 1);
            Clarfb("R", "N", "F", "C", m, n - i + 1, ib, &v[(i - 1) + (i - 1) * ldv], ldv, &t[(i - 1) * ldt], ldt, &c[(i - 1) * ldc], ldc, work, ldwork);
        }
        //
    } else if (left && notran) {
        //
        kf = ((k - 1) / nb) * nb + 1;
        for (i = kf; i >= 1; i = i - nb) {
            ib = min(nb, k - i + 1);
            Clarfb("L", "N", "F", "C", m - i + 1, n, ib, &v[(i - 1) + (i - 1) * ldv], ldv, &t[(i - 1) * ldt], ldt, &c[(i - 1)], ldc, work, ldwork);
        }
        //
    } else if (right && tran) {
        //
        kf = ((k - 1) / nb) * nb + 1;
        for (i = kf; i >= 1; i = i - nb) {
            ib = min(nb, k - i + 1);
            Clarfb("R", "C", "F", "C", m, n - i + 1, ib, &v[(i - 1) + (i - 1) * ldv], ldv, &t[(i - 1) * ldt], ldt, &c[(i - 1) * ldc], ldc, work, ldwork);
        }
        //
    }
    //
    //     End of Cgemqrt
    //
}
