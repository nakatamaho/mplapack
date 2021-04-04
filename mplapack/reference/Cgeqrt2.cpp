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

void Cgeqrt2(INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, COMPLEX *t, INTEGER const &ldt, INTEGER &info) {
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
    } else if (ldt < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Cgeqrt2", -info);
        return;
    }
    //
    INTEGER k = min(m, n);
    //
    INTEGER i = 0;
    COMPLEX aii = 0.0;
    const COMPLEX one = (1.00, 0.00);
    const COMPLEX zero = (0.00, 0.00);
    COMPLEX alpha = 0.0;
    for (i = 1; i <= k; i = i + 1) {
        //
        //        Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)
        //
        Clarfg(m - i + 1, a[(i - 1) + (i - 1) * lda], a[((min(i + 1) - 1) + (m)-1) * lda], 1, t[(i - 1)]);
        if (i < n) {
            //
            //           Apply H(i) to A(I:M,I+1:N) from the left
            //
            aii = a[(i - 1) + (i - 1) * lda];
            a[(i - 1) + (i - 1) * lda] = one;
            //
            //           W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]
            //
            Cgemv("C", m - i + 1, n - i, one, a[(i - 1) + ((i + 1) - 1) * lda], lda, a[(i - 1) + (i - 1) * lda], 1, zero, t[(n - 1) * ldt], 1);
            //
            //           A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H
            //
            alpha = -conjg[t[(i - 1)] - 1];
            Cgerc(m - i + 1, n - i, alpha, a[(i - 1) + (i - 1) * lda], 1, t[(n - 1) * ldt], 1, a[(i - 1) + ((i + 1) - 1) * lda], lda);
            a[(i - 1) + (i - 1) * lda] = aii;
        }
    }
    //
    for (i = 2; i <= n; i = i + 1) {
        aii = a[(i - 1) + (i - 1) * lda];
        a[(i - 1) + (i - 1) * lda] = one;
        //
        //        T(1:I-1,I) := alpha * A(I:M,1:I-1)**H * A(I:M,I)
        //
        alpha = -t[(i - 1)];
        Cgemv("C", m - i + 1, i - 1, alpha, a[(i - 1)], lda, a[(i - 1) + (i - 1) * lda], 1, zero, t[(i - 1) * ldt], 1);
        a[(i - 1) + (i - 1) * lda] = aii;
        //
        //        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
        //
        Ctrmv("U", "N", "N", i - 1, t, ldt, t[(i - 1) * ldt], 1);
        //
        //           T(I,I) = tau(I)
        //
        t[(i - 1) + (i - 1) * ldt] = t[(i - 1)];
        t[(i - 1)] = zero;
    }
    //
    //     End of Cgeqrt2
    //
}
