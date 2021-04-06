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

void Cunbdb5(INTEGER const m1, INTEGER const m2, INTEGER const n, COMPLEX *x1, INTEGER const incx1, COMPLEX *x2, INTEGER const incx2, COMPLEX *q1, INTEGER const ldq1, COMPLEX *q2, INTEGER const ldq2, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. External Functions ..
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
        Mxerbla("Cunbdb5", -info);
        return;
    }
    //
    //     Project X onto the orthogonal complement of Q
    //
    INTEGER childinfo = 0;
    Cunbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, childinfo);
    //
    //     If the projection is nonzero, then return
    //
    const COMPLEX zero = (0.0, 0.0);
    if (RCnrm2(m1, x1, incx1) != zero || RCnrm2(m2, x2, incx2) != zero) {
        return;
    }
    //
    //     Project each standard basis vector e_1,...,e_M1 in turn, stopping
    //     when a nonzero projection is found
    //
    INTEGER i = 0;
    INTEGER j = 0;
    const COMPLEX one = (1.0, 0.0);
    for (i = 1; i <= m1; i = i + 1) {
        for (j = 1; j <= m1; j = j + 1) {
            x1[j - 1] = zero;
        }
        x1[i - 1] = one;
        for (j = 1; j <= m2; j = j + 1) {
            x2[j - 1] = zero;
        }
        Cunbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, childinfo);
        if (RCnrm2(m1, x1, incx1) != zero || RCnrm2(m2, x2, incx2) != zero) {
            return;
        }
    }
    //
    //     Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
    //     stopping when a nonzero projection is found
    //
    for (i = 1; i <= m2; i = i + 1) {
        for (j = 1; j <= m1; j = j + 1) {
            x1[j - 1] = zero;
        }
        for (j = 1; j <= m2; j = j + 1) {
            x2[j - 1] = zero;
        }
        x2[i - 1] = one;
        Cunbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, childinfo);
        if (RCnrm2(m1, x1, incx1) != zero || RCnrm2(m2, x2, incx2) != zero) {
            return;
        }
    }
    //
    //     End of Cunbdb5
    //
}
