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

void Cungl2(INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER &info) {
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
        Mxerbla("Cungl2", -info);
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
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    if (k < m) {
        //
        //        Initialise rows k+1:m to rows of the unit matrix
        //
        for (j = 1; j <= n; j = j + 1) {
            for (l = k + 1; l <= m; l = l + 1) {
                a[(l - 1) + (j - 1) * lda] = zero;
            }
            if (j > k && j <= m) {
                a[(j - 1) + (j - 1) * lda] = one;
            }
        }
    }
    //
    INTEGER i = 0;
    for (i = k; i >= 1; i = i - 1) {
        //
        //        Apply H(i)**H to A(i:m,i:n) from the right
        //
        if (i < n) {
            Clacgv(n - i, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
            if (i < m) {
                a[(i - 1) + (i - 1) * lda] = one;
                Clarf("Right", m - i, n - i + 1, &a[(i - 1) + (i - 1) * lda], lda, conj(tau[i - 1]), &a[((i + 1) - 1) + (i - 1) * lda], lda, work);
            }
            Cscal(n - i, -tau[i - 1], &a[(i - 1) + ((i + 1) - 1) * lda], lda);
            Clacgv(n - i, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
        }
        a[(i - 1) + (i - 1) * lda] = one - conj(tau[i - 1]);
        //
        //        Set A(i,1:i-1) to zero
        //
        for (l = 1; l <= i - 1; l = l + 1) {
            a[(i - 1) + (l - 1) * lda] = zero;
        }
    }
    //
    //     End of Cungl2
    //
}
