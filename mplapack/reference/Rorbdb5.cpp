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

// Derived from LAPACK routine DORBDB5.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Rorbdb5(INTEGER const m1, INTEGER const m2, INTEGER const n, REAL *x1, INTEGER const incx1, REAL *x2, INTEGER const incx2, REAL *q1, INTEGER const ldq1, REAL *q2, INTEGER const ldq2, REAL *work, INTEGER const lwork, INTEGER &info) {
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
        Mxerbla("Rorbdb5", -info);
        return;
    }
    //
    REAL eps = Rlamch("Precision");
    //
    // Project X onto the orthogonal complement of Q if X is nonzero
    //
    const REAL realzero = 0.0;
    REAL scl = realzero;
    REAL ssq = realzero;
    Rlassq(m1, x1, incx1, scl, ssq);
    Rlassq(m2, x2, incx2, scl, ssq);
    REAL norm = scl * sqrt(ssq);
    //
    const REAL one = 1.0;
    INTEGER childinfo = 0;
    if (norm > n * eps) {
        // Scale vector to unit norm to avoid problems in the caller code.
        // Computing the reciprocal is undesirable but
        // * xLASCL cannot be used because of the vector increments and
        // * the round-off error has a negligible impact on
        // orthogonalization.
        Rscal(m1, one / norm, x1, incx1);
        Rscal(m2, one / norm, x2, incx2);
        Rorbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, childinfo);
        //
        // If the projection is nonzero, then return
        //
        if (Rnrm2(m1, x1, incx1) != realzero || Rnrm2(m2, x2, incx2) != realzero) {
            return;
        }
    }
    //
    // Project each standard basis vector e_1,...,e_M1 in turn, stopping
    // when a nonzero projection is found
    //
    INTEGER i = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    for (i = 1; i <= m1; i = i + 1) {
        for (j = 1; j <= m1; j = j + 1) {
            x1[j - 1] = zero;
        }
        x1[i - 1] = one;
        for (j = 1; j <= m2; j = j + 1) {
            x2[j - 1] = zero;
        }
        Rorbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, childinfo);
        if (Rnrm2(m1, x1, incx1) != realzero || Rnrm2(m2, x2, incx2) != realzero) {
            return;
        }
    }
    //
    // Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
    // stopping when a nonzero projection is found
    //
    for (i = 1; i <= m2; i = i + 1) {
        for (j = 1; j <= m1; j = j + 1) {
            x1[j - 1] = zero;
        }
        for (j = 1; j <= m2; j = j + 1) {
            x2[j - 1] = zero;
        }
        x2[i - 1] = one;
        Rorbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, childinfo);
        if (Rnrm2(m1, x1, incx1) != realzero || Rnrm2(m2, x2, incx2) != realzero) {
            return;
        }
    }
    //
    // End of Rorbdb5
    //
}
