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

void Cgerq2(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER &info) {
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
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Cgerq2", -info);
        return;
    }
    //
    INTEGER k = min(m, n);
    //
    INTEGER i = 0;
    COMPLEX alpha = 0.0;
    const COMPLEX one = (1.0, 0.0);
    for (i = k; i >= 1; i = i - 1) {
        //
        //        Generate elementary reflector H(i) to annihilate
        //        A(m-k+i,1:n-k+i-1)
        //
        Clacgv(n - k + i, &a[((m - k + i) - 1)], lda);
        alpha = a[((m - k + i) - 1) + ((n - k + i) - 1) * lda];
        Clarfg(n - k + i, alpha, &a[((m - k + i) - 1)], lda, tau[i - 1]);
        //
        //        Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right
        //
        a[((m - k + i) - 1) + ((n - k + i) - 1) * lda] = one;
        Clarf("Right", m - k + i - 1, n - k + i, &a[((m - k + i) - 1)], lda, tau[i - 1], a, lda, work);
        a[((m - k + i) - 1) + ((n - k + i) - 1) * lda] = alpha;
        Clacgv(n - k + i - 1, &a[((m - k + i) - 1)], lda);
    }
    //
    //     End of Cgerq2
    //
}
