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

void Rggsvp3(const char *jobu, const char *jobv, const char *jobq, INTEGER const &m, INTEGER const &p, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL const &tola, REAL const &tolb, INTEGER &k, INTEGER &l, REAL *u, INTEGER const &ldu, REAL *v, INTEGER const &ldv, REAL *q, INTEGER const &ldq, arr_ref<INTEGER> iwork, REAL *tau, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    bool wantu = Mlsame(jobu, "U");
    bool wantv = Mlsame(jobv, "V");
    bool wantq = Mlsame(jobq, "Q");
    bool forwrd = true;
    bool lquery = (lwork == -1);
    INTEGER lwkopt = 1;
    //
    //     Test the input arguments
    //
    info = 0;
    if (!(wantu || Mlsame(jobu, "N"))) {
        info = -1;
    } else if (!(wantv || Mlsame(jobv, "N"))) {
        info = -2;
    } else if (!(wantq || Mlsame(jobq, "N"))) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (p < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (lda < max((INTEGER)1, m)) {
        info = -8;
    } else if (ldb < max((INTEGER)1, p)) {
        info = -10;
    } else if (ldu < 1 || (wantu && ldu < m)) {
        info = -16;
    } else if (ldv < 1 || (wantv && ldv < p)) {
        info = -18;
    } else if (ldq < 1 || (wantq && ldq < n)) {
        info = -20;
    } else if (lwork < 1 && !lquery) {
        info = -24;
    }
    //
    //     Compute workspace
    //
    if (info == 0) {
        Rgeqp3(p, n, b, ldb, iwork, tau, work, -1, info);
        lwkopt = INTEGER(work[1 - 1]);
        if (wantv) {
            lwkopt = max(lwkopt, p);
        }
        lwkopt = max(lwkopt, min(n, p));
        lwkopt = max(lwkopt, m);
        if (wantq) {
            lwkopt = max(lwkopt, n);
        }
        Rgeqp3(m, n, a, lda, iwork, tau, work, -1, info);
        lwkopt = max(lwkopt, INTEGER(work[1 - 1]));
        lwkopt = max((INTEGER)1, lwkopt);
        work[1 - 1] = lwkopt.real();
    }
    //
    if (info != 0) {
        Mxerbla("Rggsvp3", -info);
        return;
    }
    if (lquery) {
        return;
    }
    //
    //     QR with column pivoting of B: B*P = V*( S11 S12 )
    //                                           (  0   0  )
    //
    INTEGER i = 0;
    for (i = 1; i <= n; i = i + 1) {
        iwork[i - 1] = 0;
    }
    Rgeqp3(p, n, b, ldb, iwork, tau, work, lwork, info);
    //
    //     Update A := A*P
    //
    Rlapmt(forwrd, m, n, a, lda, iwork);
    //
    //     Determine the effective rank of matrix B.
    //
    l = 0;
    for (i = 1; i <= min(p, n); i = i + 1) {
        if (abs(b[(i - 1) + (i - 1) * ldb]) > tolb) {
            l++;
        }
    }
    //
    const REAL zero = 0.0;
    if (wantv) {
        //
        //        Copy the details of V, and form V.
        //
        Rlaset("Full", p, p, zero, zero, v, ldv);
        if (p > 1) {
            Rlacpy("Lower", p - 1, n, b[(2 - 1)], ldb, v[(2 - 1)], ldv);
        }
        Rorg2r(p, p, min(p, n), v, ldv, tau, work, info);
    }
    //
    //     Clean up B
    //
    INTEGER j = 0;
    for (j = 1; j <= l - 1; j = j + 1) {
        for (i = j + 1; i <= l; i = i + 1) {
            b[(i - 1) + (j - 1) * ldb] = zero;
        }
    }
    if (p > l) {
        Rlaset("Full", p - l, n, zero, zero, b[((l + 1) - 1)], ldb);
    }
    //
    const REAL one = 1.0;
    if (wantq) {
        //
        //        Set Q = I and Update Q := Q*P
        //
        Rlaset("Full", n, n, zero, one, q, ldq);
        Rlapmt(forwrd, n, n, q, ldq, iwork);
    }
    //
    if (p >= l && n != l) {
        //
        //        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z
        //
        Rgerq2(l, n, b, ldb, tau, work, info);
        //
        //        Update A := A*Z**T
        //
        Rormr2("Right", "Transpose", m, n, l, b, ldb, tau, a, lda, work, info);
        //
        if (wantq) {
            //
            //           Update Q := Q*Z**T
            //
            Rormr2("Right", "Transpose", n, n, l, b, ldb, tau, q, ldq, work, info);
        }
        //
        //        Clean up B
        //
        Rlaset("Full", l, n - l, zero, zero, b, ldb);
        for (j = n - l + 1; j <= n; j = j + 1) {
            for (i = j - n + l + 1; i <= l; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = zero;
            }
        }
        //
    }
    //
    //     Let              N-L     L
    //                A = ( A11    A12 ) M,
    //
    //     then the following does the complete QR decomposition of A11:
    //
    //              A11 = U*(  0  T12 )*P1**T
    //                      (  0   0  )
    //
    for (i = 1; i <= n - l; i = i + 1) {
        iwork[i - 1] = 0;
    }
    Rgeqp3(m, n - l, a, lda, iwork, tau, work, lwork, info);
    //
    //     Determine the effective rank of A11
    //
    k = 0;
    for (i = 1; i <= min(m, n - l); i = i + 1) {
        if (abs(a[(i - 1) + (i - 1) * lda]) > tola) {
            k++;
        }
    }
    //
    //     Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N )
    //
    Rorm2r("Left", "Transpose", m, l, min(m, n - l), a, lda, tau, a[((n - l + 1) - 1) * lda], lda, work, info);
    //
    if (wantu) {
        //
        //        Copy the details of U, and form U
        //
        Rlaset("Full", m, m, zero, zero, u, ldu);
        if (m > 1) {
            Rlacpy("Lower", m - 1, n - l, a[(2 - 1)], lda, u[(2 - 1)], ldu);
        }
        Rorg2r(m, m, min(m, n - l), u, ldu, tau, work, info);
    }
    //
    if (wantq) {
        //
        //        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
        //
        Rlapmt(forwrd, n, n - l, q, ldq, iwork);
    }
    //
    //     Clean up A: set the strictly lower triangular part of
    //     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.
    //
    for (j = 1; j <= k - 1; j = j + 1) {
        for (i = j + 1; i <= k; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
    }
    if (m > k) {
        Rlaset("Full", m - k, n - l, zero, zero, a[((k + 1) - 1)], lda);
    }
    //
    if (n - l > k) {
        //
        //        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
        //
        Rgerq2(k, n - l, a, lda, tau, work, info);
        //
        if (wantq) {
            //
            //           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T
            //
            Rormr2("Right", "Transpose", n, n - l, k, a, lda, tau, q, ldq, work, info);
        }
        //
        //        Clean up A
        //
        Rlaset("Full", k, n - l - k, zero, zero, a, lda);
        for (j = n - l - k + 1; j <= n - l; j = j + 1) {
            for (i = j - n + l + k + 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = zero;
            }
        }
        //
    }
    //
    if (m > k) {
        //
        //        QR factorization of A( K+1:M,N-L+1:N )
        //
        Rgeqr2(m - k, l, a[((k + 1) - 1) + ((n - l + 1) - 1) * lda], lda, tau, work, info);
        //
        if (wantu) {
            //
            //           Update U(:,K+1:M) := U(:,K+1:M)*U1
            //
            Rorm2r("Right", "No transpose", m, m - k, min(m - k, l), a[((k + 1) - 1) + ((n - l + 1) - 1) * lda], lda, tau, u[((k + 1) - 1) * ldu], ldu, work, info);
        }
        //
        //        Clean up
        //
        for (j = n - l + 1; j <= n; j = j + 1) {
            for (i = j - n + k + l + 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = zero;
            }
        }
        //
    }
    //
    work[1 - 1] = lwkopt.real();
    //
    //     End of Rggsvp3
    //
}
