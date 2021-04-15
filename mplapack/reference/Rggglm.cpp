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

void Rggglm(INTEGER const n, INTEGER const m, INTEGER const p, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *d, REAL *x, REAL *y, REAL *work, INTEGER const lwork, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  ===================================================================
    //
    //     .. Parameters ..
    //     ..
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
    INTEGER np = min(n, p);
    bool lquery = (lwork == -1);
    if (n < 0) {
        info = -1;
    } else if (m < 0 || m > n) {
        info = -2;
    } else if (p < 0 || p < n - m) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    }
    //
    //     Calculate workspace
    //
    INTEGER lwkmin = 0;
    INTEGER lwkopt = 0;
    INTEGER nb1 = 0;
    INTEGER nb2 = 0;
    INTEGER nb3 = 0;
    INTEGER nb4 = 0;
    INTEGER nb = 0;
    if (info == 0) {
        if (n == 0) {
            lwkmin = 1;
            lwkopt = 1;
        } else {
            nb1 = iMlaenv(1, "Rgeqrf", " ", n, m, -1, -1);
            nb2 = iMlaenv(1, "RgerQF", " ", n, m, -1, -1);
            nb3 = iMlaenv(1, "Rormqr", " ", n, m, p, -1);
            nb4 = iMlaenv(1, "Rormrq", " ", n, m, p, -1);
            nb = max(nb1, nb2, nb3, nb4);
            lwkmin = m + n + p;
            lwkopt = m + np + max(n, p) * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < lwkmin && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rggglm", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    if (n == 0) {
        for (i = 1; i <= m; i = i + 1) {
            x[i - 1] = zero;
        }
        for (i = 1; i <= p; i = i + 1) {
            y[i - 1] = zero;
        }
        return;
    }
    //
    //     Compute the GQR factorization of matrices A and B:
    //
    //          Q**T*A = ( R11 ) M,    Q**T*B*Z**T = ( T11   T12 ) M
    //                   (  0  ) N-M                 (  0    T22 ) N-M
    //                      M                         M+P-N  N-M
    //
    //     where R11 and T22 are upper triangular, and Q and Z are
    //     orthogonal.
    //
    Rggqrf(n, m, p, a, lda, work, b, ldb, &work[(m + 1) - 1], &work[(m + np + 1) - 1], lwork - m - np, info);
    INTEGER lopt = work[(m + np + 1) - 1];
    //
    //     Update left-hand-side vector d = Q**T*d = ( d1 ) M
    //                                               ( d2 ) N-M
    //
    Rormqr("Left", "Transpose", n, 1, m, a, lda, work, d, max((INTEGER)1, n), &work[(m + np + 1) - 1], lwork - m - np, info);
    lopt = max(lopt, int(work[(m + np + 1) - 1]));
    //
    //     Solve T22*y2 = d2 for y2
    //
    if (n > m) {
        Rtrtrs("Upper", "No transpose", "Non unit", n - m, 1, &b[((m + 1) - 1) + ((m + p - n + 1) - 1) * ldb], ldb, &d[(m + 1) - 1], n - m, info);
        //
        if (info > 0) {
            info = 1;
            return;
        }
        //
        Rcopy(n - m, &d[(m + 1) - 1], 1, y[(m + p - n + 1) - 1], 1);
    }
    //
    //     Set y1 = 0
    //
    for (i = 1; i <= m + p - n; i = i + 1) {
        y[i - 1] = zero;
    }
    //
    //     Update d1 = d1 - T12*y2
    //
    const REAL one = 1.0;
    Rgemv("No transpose", m, n - m, -one, &b[((m + p - n + 1) - 1) * ldb], ldb, y[(m + p - n + 1) - 1], 1, one, d, 1);
    //
    //     Solve triangular system: R11*x = d1
    //
    if (m > 0) {
        Rtrtrs("Upper", "No Transpose", "Non unit", m, 1, a, lda, d, m, info);
        //
        if (info > 0) {
            info = 2;
            return;
        }
        //
        //        Copy D to X
        //
        Rcopy(m, d, 1, x, 1);
    }
    //
    //     Backward transformation y = Z**T *y
    //
    Rormrq("Left", "Transpose", p, 1, np, &b[(max((n - p + 1)) - 1) * ldb], ldb, &work[(m + 1) - 1], y, max((INTEGER)1, p), &work[(m + np + 1) - 1], lwork - m - np, info);
    work[1 - 1] = m + np + max(lopt, int(work[(m + np + 1) - 1]));
    //
    //     End of Rggglm
    //
}
