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

void Rorbdb6(INTEGER const m1, INTEGER const m2, INTEGER const n, REAL *x1, INTEGER const incx1, REAL *x2, INTEGER const incx2, REAL *q1, INTEGER const ldq1, REAL *q2, INTEGER const ldq2, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Function ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test input arguments
    //
    info = 0;
    if (m1 < 0) {
        info = -1;
    } else if (m2 < 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (incx1 < 1) {
        info = -5;
    } else if (incx2 < 1) {
        info = -7;
    } else if (ldq1 < max((INTEGER)1, m1)) {
        info = -9;
    } else if (ldq2 < max((INTEGER)1, m2)) {
        info = -11;
    } else if (lwork < n) {
        info = -13;
    }
    //
    if (info != 0) {
        Mxerbla("Rorbdb6", -info);
        return;
    }
    //
    //     First, project X onto the orthogonal complement of Q's column
    //     space
    //
    const REAL realzero = 0.0;
    REAL scl1 = realzero;
    const REAL realone = 1.0;
    REAL ssq1 = realone;
    Rlassq(m1, x1, incx1, scl1, ssq1);
    REAL scl2 = realzero;
    REAL ssq2 = realone;
    Rlassq(m2, x2, incx2, scl2, ssq2);
    REAL normsq1 = pow2(scl1) * ssq1 + pow2(scl2) * ssq2;
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (m1 == 0) {
        for (i = 1; i <= n; i = i + 1) {
            work[i - 1] = zero;
        }
    } else {
        Rgemv("C", m1, n, one, q1, ldq1, x1, incx1, zero, work, 1);
    }
    //
    Rgemv("C", m2, n, one, q2, ldq2, x2, incx2, one, work, 1);
    //
    const REAL negone = -1.0;
    Rgemv("N", m1, n, negone, q1, ldq1, work, 1, one, x1, incx1);
    Rgemv("N", m2, n, negone, q2, ldq2, work, 1, one, x2, incx2);
    //
    scl1 = realzero;
    ssq1 = realone;
    Rlassq(m1, x1, incx1, scl1, ssq1);
    scl2 = realzero;
    ssq2 = realone;
    Rlassq(m2, x2, incx2, scl2, ssq2);
    REAL normsq2 = pow2(scl1) * ssq1 + pow2(scl2) * ssq2;
    //
    //     If projection is sufficiently large in norm, then stop.
    //     If projection is zero, then stop.
    //     Otherwise, project again.
    //
    const REAL alphasq = 0.01e0;
    if (normsq2 >= alphasq * normsq1) {
        return;
    }
    //
    if (normsq2 == zero) {
        return;
    }
    //
    normsq1 = normsq2;
    //
    for (i = 1; i <= n; i = i + 1) {
        work[i - 1] = zero;
    }
    //
    if (m1 == 0) {
        for (i = 1; i <= n; i = i + 1) {
            work[i - 1] = zero;
        }
    } else {
        Rgemv("C", m1, n, one, q1, ldq1, x1, incx1, zero, work, 1);
    }
    //
    Rgemv("C", m2, n, one, q2, ldq2, x2, incx2, one, work, 1);
    //
    Rgemv("N", m1, n, negone, q1, ldq1, work, 1, one, x1, incx1);
    Rgemv("N", m2, n, negone, q2, ldq2, work, 1, one, x2, incx2);
    //
    scl1 = realzero;
    ssq1 = realone;
    Rlassq(m1, x1, incx1, scl1, ssq1);
    scl2 = realzero;
    ssq2 = realone;
    Rlassq(m1, x1, incx1, scl1, ssq1);
    normsq2 = pow2(scl1) * ssq1 + pow2(scl2) * ssq2;
    //
    //     If second projection is sufficiently large in norm, then do
    //     nothing more. Alternatively, if it shrunk significantly, then
    //     truncate it to zero.
    //
    if (normsq2 < alphasq * normsq1) {
        for (i = 1; i <= m1; i = i + 1) {
            x1[i - 1] = zero;
        }
        for (i = 1; i <= m2; i = i + 1) {
            x2[i - 1] = zero;
        }
    }
    //
    //     End of Rorbdb6
    //
}
