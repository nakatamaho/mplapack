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

void Rorgr2(INTEGER const m, INTEGER const n, INTEGER const k, REAL *a, INTEGER const lda, REAL *tau, REAL *work, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < m) {
        info = -2;
    } else if (k < 0 || k > m) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Rorgr2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m <= 0) {
        return;
    }
    //
    INTEGER j = 0;
    INTEGER l = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (k < m) {
        //
        //        Initialise rows 1:m-k to rows of the unit matrix
        //
        for (j = 1; j <= n; j = j + 1) {
            for (l = 1; l <= m - k; l = l + 1) {
                a[(l - 1) + (j - 1) * lda] = zero;
            }
            if (j > n - m && j <= n - k) {
                a[((m - n + j) - 1) + (j - 1) * lda] = one;
            }
        }
    }
    //
    INTEGER i = 0;
    INTEGER ii = 0;
    for (i = 1; i <= k; i = i + 1) {
        ii = m - k + i;
        //
        //        Apply H(i) to A(1:m-k+i,1:n-k+i) from the right
        //
        a[(ii - 1) + ((n - m + ii) - 1) * lda] = one;
        Rlarf("Right", ii - 1, n - m + ii, &a[(ii - 1)], lda, tau[i - 1], a, lda, work);
        Rscal(n - m + ii - 1, -tau[i - 1], &a[(ii - 1)], lda);
        a[(ii - 1) + ((n - m + ii) - 1) * lda] = one - tau[i - 1];
        //
        //        Set A(m-k+i,n-k+i+1:n) to zero
        //
        for (l = n - m + ii + 1; l <= n; l = l + 1) {
            a[(ii - 1) + (l - 1) * lda] = zero;
        }
    }
    //
    //     End of Rorgr2
    //
}
