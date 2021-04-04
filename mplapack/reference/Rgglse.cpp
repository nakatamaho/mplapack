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

void Rgglse(INTEGER const &m, INTEGER const &n, INTEGER const &p, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *c, REAL *d, REAL *x, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    INTEGER mn = min(m, n);
    bool lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (p < 0 || p > n || p < n - m) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, p)) {
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
            nb1 = iMlaenv[("Rgeqrf" - 1) * ldiMlaenv];
            nb2 = iMlaenv[("Rgerqf" - 1) * ldiMlaenv];
            nb3 = iMlaenv[("Rormqr" - 1) * ldiMlaenv];
            nb4 = iMlaenv[("Rormrq" - 1) * ldiMlaenv];
            nb = max(nb1, nb2, nb3, nb4);
            lwkmin = m + n + p;
            lwkopt = p + mn + max(m, n) * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < lwkmin && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgglse", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Compute the GRQ factorization of matrices B and A:
    //
    //            B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P
    //                        N-P  P                     (  0  R22 ) M+P-N
    //                                                      N-P  P
    //
    //     where T12 and R11 are upper triangular, and Q and Z are
    //     orthogonal.
    //
    Rggrqf(p, m, n, b, ldb, work, a, lda, work[(p + 1) - 1], work[(p + mn + 1) - 1], lwork - p - mn, info);
    INTEGER lopt = work[(p + mn + 1) - 1];
    //
    //     Update c = Z**T *c = ( c1 ) N-P
    //                          ( c2 ) M+P-N
    //
    Rormqr("Left", "Transpose", m, 1, mn, a, lda, work[(p + 1) - 1], c, max((INTEGER)1, m), work[(p + mn + 1) - 1], lwork - p - mn, info);
    lopt = max(lopt, INTEGER(work[(p + mn + 1) - 1]));
    //
    //     Solve T12*x2 = d for x2
    //
    const REAL one = 1.0;
    if (p > 0) {
        Rtrtrs("Upper", "No transpose", "Non-unit", p, 1, b[((n - p + 1) - 1) * ldb], ldb, d, p, info);
        //
        if (info > 0) {
            info = 1;
            return;
        }
        //
        //        Put the solution in X
        //
        Rcopy(p, d, 1, x[(n - p + 1) - 1], 1);
        //
        //        Update c1
        //
        Rgemv("No transpose", n - p, p, -one, a[((n - p + 1) - 1) * lda], lda, d, 1, one, c, 1);
    }
    //
    //     Solve R11*x1 = c1 for x1
    //
    if (n > p) {
        Rtrtrs("Upper", "No transpose", "Non-unit", n - p, 1, a, lda, c, n - p, info);
        //
        if (info > 0) {
            info = 2;
            return;
        }
        //
        //        Put the solutions in X
        //
        Rcopy(n - p, c, 1, x, 1);
    }
    //
    //     Compute the residual vector:
    //
    INTEGER nr = 0;
    if (m < n) {
        nr = m + p - n;
        if (nr > 0) {
            Rgemv("No transpose", nr, n - m, -one, a[((n - p + 1) - 1) + ((m + 1) - 1) * lda], lda, d[(nr + 1) - 1], 1, one, c[(n - p + 1) - 1], 1);
        }
    } else {
        nr = p;
    }
    if (nr > 0) {
        Rtrmv("Upper", "No transpose", "Non unit", nr, a[((n - p + 1) - 1) + ((n - p + 1) - 1) * lda], lda, d, 1);
        Raxpy(nr, -one, d, 1, c[(n - p + 1) - 1], 1);
    }
    //
    //     Backward transformation x = Q**T*x
    //
    Rormrq("Left", "Transpose", n, 1, p, b, ldb, work[1 - 1], x, n, work[(p + mn + 1) - 1], lwork - p - mn, info);
    work[1 - 1] = p + mn + max(lopt, INTEGER(work[(p + mn + 1) - 1]));
    //
    //     End of Rgglse
    //
}
