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

// Derived from LAPACK routine ZUNBDB6.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Cunbdb6(INTEGER const m1, INTEGER const m2, INTEGER const n, COMPLEX *x1, INTEGER const incx1, COMPLEX *x2, INTEGER const incx2, COMPLEX *q1, INTEGER const ldq1, COMPLEX *q2, INTEGER const ldq2, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
    //
    // Test input arguments
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
        Mxerbla("Cunbdb6", -info);
        return;
    }
    //
    REAL eps = Rlamch("Precision");
    //
    // Compute the Euclidean norm of X
    //
    const REAL realzero = 0.0;
    REAL scl = realzero;
    REAL ssq = realzero;
    Classq(m1, x1, incx1, scl, ssq);
    Classq(m2, x2, incx2, scl, ssq);
    REAL norm = scl * sqrt(ssq);
    //
    // First, project X onto the orthogonal complement of Q's column
    // space
    //
    INTEGER i = 0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    if (m1 == 0) {
        for (i = 1; i <= n; i = i + 1) {
            work[i - 1] = zero;
        }
    } else {
        Cgemv("C", m1, n, one, q1, ldq1, x1, incx1, zero, work, 1);
    }
    //
    Cgemv("C", m2, n, one, q2, ldq2, x2, incx2, one, work, 1);
    //
    const COMPLEX negone = COMPLEX(-1.0, 0.0);
    Cgemv("N", m1, n, negone, q1, ldq1, work, 1, one, x1, incx1);
    Cgemv("N", m2, n, negone, q2, ldq2, work, 1, one, x2, incx2);
    //
    scl = realzero;
    ssq = realzero;
    Classq(m1, x1, incx1, scl, ssq);
    Classq(m2, x2, incx2, scl, ssq);
    REAL norm_new = scl * sqrt(ssq);
    //
    // If projection is sufficiently large in norm, then stop.
    // If projection is zero, then stop.
    // Otherwise, project again.
    //
    const REAL alpha = 0.83;
    if (norm_new >= alpha * norm) {
        return;
    }
    //
    INTEGER ix = 0;
    if (norm_new <= n * eps * norm) {
        for (ix = 1; ix <= 1 + (m1 - 1) * incx1; ix = ix + incx1) {
            x1[ix - 1] = zero;
        }
        for (ix = 1; ix <= 1 + (m2 - 1) * incx2; ix = ix + incx2) {
            x2[ix - 1] = zero;
        }
        return;
    }
    //
    norm = norm_new;
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
        Cgemv("C", m1, n, one, q1, ldq1, x1, incx1, zero, work, 1);
    }
    //
    Cgemv("C", m2, n, one, q2, ldq2, x2, incx2, one, work, 1);
    //
    Cgemv("N", m1, n, negone, q1, ldq1, work, 1, one, x1, incx1);
    Cgemv("N", m2, n, negone, q2, ldq2, work, 1, one, x2, incx2);
    //
    scl = realzero;
    ssq = realzero;
    Classq(m1, x1, incx1, scl, ssq);
    Classq(m2, x2, incx2, scl, ssq);
    norm_new = scl * sqrt(ssq);
    //
    // If second projection is sufficiently large in norm, then do
    // nothing more. Alternatively, if it shrunk significantly, then
    // truncate it to zero.
    //
    if (norm_new < alpha * norm) {
        for (ix = 1; ix <= 1 + (m1 - 1) * incx1; ix = ix + incx1) {
            x1[ix - 1] = zero;
        }
        for (ix = 1; ix <= 1 + (m2 - 1) * incx2; ix = ix + incx2) {
            x2[ix - 1] = zero;
        }
    }
    //
    // End of Cunbdb6
    //
}
