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

void Rorgbr(const char *vect, INTEGER const m, INTEGER const n, INTEGER const k, REAL *a, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     Test the input arguments
    //
    info = 0;
    bool wantq = Mlsame(vect, "Q");
    INTEGER mn = min(m, n);
    bool lquery = (lwork == -1);
    if (!wantq && !Mlsame(vect, "P")) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (n < 0 || (wantq && (n > m || n < min(m, k))) || (!wantq && (m > n || m < min(n, k)))) {
        info = -3;
    } else if (k < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (lwork < max((INTEGER)1, mn) && !lquery) {
        info = -9;
    }
    //
    INTEGER iinfo = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        work[1 - 1] = 1;
        if (wantq) {
            if (m >= k) {
                Rorgqr(m, n, k, a, lda, tau, work, -1, iinfo);
            } else {
                if (m > 1) {
                    Rorgqr(m - 1, m - 1, m - 1, &a[(2 - 1) + (2 - 1) * lda], lda, tau, work, -1, iinfo);
                }
            }
        } else {
            if (k < n) {
                Rorglq(m, n, k, a, lda, tau, work, -1, iinfo);
            } else {
                if (n > 1) {
                    Rorglq(n - 1, n - 1, n - 1, &a[(2 - 1) + (2 - 1) * lda], lda, tau, work, -1, iinfo);
                }
            }
        }
        lwkopt = work[1 - 1];
        lwkopt = max(lwkopt, mn);
    }
    //
    if (info != 0) {
        Mxerbla("Rorgbr", -info);
        return;
    } else if (lquery) {
        work[1 - 1] = lwkopt;
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    INTEGER j = 0;
    const REAL zero = 0.0;
    INTEGER i = 0;
    const REAL one = 1.0;
    if (wantq) {
        //
        //        Form Q, determined by a call to Rgebrd to reduce an m-by-k
        //        matrix
        //
        if (m >= k) {
            //
            //           If m >= k, assume m >= n >= k
            //
            Rorgqr(m, n, k, a, lda, tau, work, lwork, iinfo);
            //
        } else {
            //
            //           If m < k, assume m = n
            //
            //           Shift the vectors which define the elementary reflectors one
            //           column to the right, and set the first row and column of Q
            //           to those of the unit matrix
            //
            for (j = m; j >= 2; j = j - 1) {
                a[(j - 1) * lda] = zero;
                for (i = j + 1; i <= m; i = i + 1) {
                    a[(i - 1) + (j - 1) * lda] = a[(i - 1) + ((j - 1) - 1) * lda];
                }
            }
            a[(1 - 1)] = one;
            for (i = 2; i <= m; i = i + 1) {
                a[(i - 1)] = zero;
            }
            if (m > 1) {
                //
                //              Form Q(2:m,2:m)
                //
                Rorgqr(m - 1, m - 1, m - 1, &a[(2 - 1) + (2 - 1) * lda], lda, tau, work, lwork, iinfo);
            }
        }
    } else {
        //
        //        Form P**T, determined by a call to Rgebrd to reduce a k-by-n
        //        matrix
        //
        if (k < n) {
            //
            //           If k < n, assume k <= m <= n
            //
            Rorglq(m, n, k, a, lda, tau, work, lwork, iinfo);
            //
        } else {
            //
            //           If k >= n, assume m = n
            //
            //           Shift the vectors which define the elementary reflectors one
            //           row downward, and set the first row and column of P**T to
            //           those of the unit matrix
            //
            a[(1 - 1)] = one;
            for (i = 2; i <= n; i = i + 1) {
                a[(i - 1)] = zero;
            }
            for (j = 2; j <= n; j = j + 1) {
                for (i = j - 1; i >= 2; i = i - 1) {
                    a[(i - 1) + (j - 1) * lda] = a[((i - 1) - 1) + (j - 1) * lda];
                }
                a[(j - 1) * lda] = zero;
            }
            if (n > 1) {
                //
                //              Form P**T(2:n,2:n)
                //
                Rorglq(n - 1, n - 1, n - 1, &a[(2 - 1) + (2 - 1) * lda], lda, tau, work, lwork, iinfo);
            }
        }
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Rorgbr
    //
}
