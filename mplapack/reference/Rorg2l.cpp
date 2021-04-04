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

void Rorg2l(INTEGER const &m, INTEGER const &n, INTEGER const &k, REAL *a, INTEGER const &lda, REAL *tau, REAL *work, INTEGER &info) {
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
    } else if (n < 0 || n > m) {
        info = -2;
    } else if (k < 0 || k > n) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Rorg2l", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     Initialise columns 1:n-k to columns of the unit matrix
    //
    INTEGER j = 0;
    INTEGER l = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    for (j = 1; j <= n - k; j = j + 1) {
        for (l = 1; l <= m; l = l + 1) {
            a[(l - 1) + (j - 1) * lda] = zero;
        }
        a[((m - n + j) - 1) + (j - 1) * lda] = one;
    }
    //
    INTEGER i = 0;
    INTEGER ii = 0;
    for (i = 1; i <= k; i = i + 1) {
        ii = n - k + i;
        //
        //        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
        //
        a[((m - n + ii) - 1) + (ii - 1) * lda] = one;
        Rlarf("Left", m - n + ii, ii - 1, a[(ii - 1) * lda], 1, tau[i - 1], a, lda, work);
        Rscal(m - n + ii - 1, -tau[i - 1], a[(ii - 1) * lda], 1);
        a[((m - n + ii) - 1) + (ii - 1) * lda] = one - tau[i - 1];
        //
        //        Set A(m-k+i+1:m,n-k+i) to zero
        //
        for (l = m - n + ii + 1; l <= m; l = l + 1) {
            a[(l - 1) + (ii - 1) * lda] = zero;
        }
    }
    //
    //     End of Rorg2l
    //
}
