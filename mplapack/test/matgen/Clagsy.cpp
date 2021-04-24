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

void Clagsy(INTEGER const n, INTEGER const k, REAL *d, COMPLEX *a, INTEGER const lda, INTEGER *iseed, COMPLEX *work, INTEGER &info) {
    //
    //  -- LAPACK auxiliary routine --
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
    //     Test the input arguments
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (k < 0 || k > n - 1) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    }
    if (info < 0) {
        Mxerbla("Clagsy", -info);
        return;
    }
    //
    //     initialize lower triangle of A to diagonal matrix
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    for (j = 1; j <= n; j = j + 1) {
        for (i = j + 1; i <= n; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        a[(i - 1) + (i - 1) * lda] = d[i - 1];
    }
    //
    //     Generate lower triangle of symmetric matrix
    //
    REAL wn = 0.0;
    COMPLEX wa = 0.0;
    COMPLEX tau = 0.0;
    COMPLEX wb = 0.0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    const COMPLEX half = COMPLEX(0.5e+0, 0.0);
    COMPLEX alpha = 0.0;
    INTEGER jj = 0;
    INTEGER ii = 0;
    for (i = n - 1; i >= 1; i = i - 1) {
        //
        //        generate random reflection
        //
        zlarnv(3, iseed, n - i + 1, work);
        wn = RCnrm2(n - i + 1, work, 1);
        wa = (wn / abs(work[1 - 1])) * work[1 - 1];
        if (wn == zero) {
            tau = zero;
        } else {
            wb = work[1 - 1] + wa;
            Cscal(n - i, one / wb, &work[2 - 1], 1);
            work[1 - 1] = one;
            tau = (wb / wa).real();
        }
        //
        //        apply random reflection to A(i:n,i:n) from the left
        //        and the right
        //
        //        compute  y := tau * A * conj(u)
        //
        zlacgv(n - i + 1, work, 1);
        zsymv("Lower", n - i + 1, tau, &a[(i - 1) + (i - 1) * lda], lda, work, 1, zero, &work[(n + 1) - 1], 1);
        zlacgv(n - i + 1, work, 1);
        //
        //        compute  v := y - 1/2 * tau * ( u, y ) * u
        //
        alpha = -half * tau * Cdotc(n - i + 1, work, 1, &work[(n + 1) - 1], 1);
        Caxpy(n - i + 1, alpha, work, 1, &work[(n + 1) - 1], 1);
        //
        //        apply the transformation as a rank-2 update to A(i:n,i:n)
        //
        //        CALL ZSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1,
        //        $               A( I, I ), LDA )
        //
        for (jj = i; jj <= n; jj = jj + 1) {
            for (ii = jj; ii <= n; ii = ii + 1) {
                a[(ii - 1) + (jj - 1) * lda] = a[(ii - 1) + (jj - 1) * lda] - work[(ii - i + 1) - 1] * work[(n + jj - i + 1) - 1] - work[(n + ii - i + 1) - 1] * work[(jj - i + 1) - 1];
            }
        }
    }
    //
    //     Reduce number of subdiagonals to K
    //
    for (i = 1; i <= n - 1 - k; i = i + 1) {
        //
        //        generate reflection to annihilate A(k+i+1:n,i)
        //
        wn = RCnrm2(n - k - i + 1, &a[((k + i) - 1) + (i - 1) * lda], 1);
        wa = (wn / abs(a[((k + i) - 1) + (i - 1) * lda])) * a[((k + i) - 1) + (i - 1) * lda];
        if (wn == zero) {
            tau = zero;
        } else {
            wb = a[((k + i) - 1) + (i - 1) * lda] + wa;
            Cscal(n - k - i, one / wb, &a[((k + i + 1) - 1) + (i - 1) * lda], 1);
            a[((k + i) - 1) + (i - 1) * lda] = one;
            tau = (wb / wa).real();
        }
        //
        //        apply reflection to A(k+i:n,i+1:k+i-1) from the left
        //
        Cgemv("Conjugate transpose", n - k - i + 1, k - 1, one, &a[((k + i) - 1) + ((i + 1) - 1) * lda], lda, &a[((k + i) - 1) + (i - 1) * lda], 1, zero, work, 1);
        Cgerc(n - k - i + 1, k - 1, -tau, &a[((k + i) - 1) + (i - 1) * lda], 1, work, 1, &a[((k + i) - 1) + ((i + 1) - 1) * lda], lda);
        //
        //        apply reflection to A(k+i:n,k+i:n) from the left and the right
        //
        //        compute  y := tau * A * conj(u)
        //
        zlacgv(n - k - i + 1, &a[((k + i) - 1) + (i - 1) * lda], 1);
        zsymv("Lower", n - k - i + 1, tau, &a[((k + i) - 1) + ((k + i) - 1) * lda], lda, &a[((k + i) - 1) + (i - 1) * lda], 1, zero, work, 1);
        zlacgv(n - k - i + 1, &a[((k + i) - 1) + (i - 1) * lda], 1);
        //
        //        compute  v := y - 1/2 * tau * ( u, y ) * u
        //
        alpha = -half * tau * Cdotc(n - k - i + 1, &a[((k + i) - 1) + (i - 1) * lda], 1, work, 1);
        Caxpy(n - k - i + 1, alpha, &a[((k + i) - 1) + (i - 1) * lda], 1, work, 1);
        //
        //        apply symmetric rank-2 update to A(k+i:n,k+i:n)
        //
        //        CALL ZSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1,
        //        $               A( K+I, K+I ), LDA )
        //
        for (jj = k + i; jj <= n; jj = jj + 1) {
            for (ii = jj; ii <= n; ii = ii + 1) {
                a[(ii - 1) + (jj - 1) * lda] = a[(ii - 1) + (jj - 1) * lda] - a[(ii - 1) + (i - 1) * lda] * work[(jj - k - i + 1) - 1] - work[(ii - k - i + 1) - 1] * a[(jj - 1) + (i - 1) * lda];
            }
        }
        //
        a[((k + i) - 1) + (i - 1) * lda] = -wa;
        for (j = k + i + 1; j <= n; j = j + 1) {
            a[(j - 1) + (i - 1) * lda] = zero;
        }
    }
    //
    //     Store full symmetric matrix
    //
    for (j = 1; j <= n; j = j + 1) {
        for (i = j + 1; i <= n; i = i + 1) {
            a[(j - 1) + (i - 1) * lda] = a[(i - 1) + (j - 1) * lda];
        }
    }
    //
    //     End of Clagsy
    //
}
