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

void Rlatrd(const char *uplo, INTEGER const &n, INTEGER const &nb, REAL *a, INTEGER const &lda, REAL *e, REAL *tau, REAL *w, INTEGER const &ldw) {
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
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    INTEGER i = 0;
    INTEGER iw = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    const REAL half = 0.5e+0;
    REAL alpha = 0.0;
    if (Mlsame(uplo, "U")) {
        //
        //        Reduce last NB columns of upper triangle
        //
        for (i = n; i >= n - nb + 1; i = i - 1) {
            iw = i - n + nb;
            if (i < n) {
                //
                //              Update A(1:i,i)
                //
                Rgemv("No transpose", i, n - i, -one, a[((i + 1) - 1) * lda], lda, w[(i - 1) + ((iw + 1) - 1) * ldw], ldw, one, a[(i - 1) * lda], 1);
                Rgemv("No transpose", i, n - i, -one, w[((iw + 1) - 1) * ldw], ldw, a[(i - 1) + ((i + 1) - 1) * lda], lda, one, a[(i - 1) * lda], 1);
            }
            if (i > 1) {
                //
                //              Generate elementary reflector H(i) to annihilate
                //              A(1:i-2,i)
                //
                Rlarfg(i - 1, a[((i - 1) - 1) + (i - 1) * lda], a[(i - 1) * lda], 1, tau[(i - 1) - 1]);
                e[(i - 1) - 1] = a[((i - 1) - 1) + (i - 1) * lda];
                a[((i - 1) - 1) + (i - 1) * lda] = one;
                //
                //              Compute W(1:i-1,i)
                //
                Rsymv("Upper", i - 1, one, a, lda, a[(i - 1) * lda], 1, zero, w[(iw - 1) * ldw], 1);
                if (i < n) {
                    Rgemv("Transpose", i - 1, n - i, one, w[((iw + 1) - 1) * ldw], ldw, a[(i - 1) * lda], 1, zero, w[((i + 1) - 1) + (iw - 1) * ldw], 1);
                    Rgemv("No transpose", i - 1, n - i, -one, a[((i + 1) - 1) * lda], lda, w[((i + 1) - 1) + (iw - 1) * ldw], 1, one, w[(iw - 1) * ldw], 1);
                    Rgemv("Transpose", i - 1, n - i, one, a[((i + 1) - 1) * lda], lda, a[(i - 1) * lda], 1, zero, w[((i + 1) - 1) + (iw - 1) * ldw], 1);
                    Rgemv("No transpose", i - 1, n - i, -one, w[((iw + 1) - 1) * ldw], ldw, w[((i + 1) - 1) + (iw - 1) * ldw], 1, one, w[(iw - 1) * ldw], 1);
                }
                Rscal(i - 1, tau[(i - 1) - 1], w[(iw - 1) * ldw], 1);
                alpha = -half * tau[(i - 1) - 1] * Rdot(i - 1, w[(iw - 1) * ldw], 1, a[(i - 1) * lda], 1);
                Raxpy(i - 1, alpha, a[(i - 1) * lda], 1, w[(iw - 1) * ldw], 1);
            }
            //
        }
    } else {
        //
        //        Reduce first NB columns of lower triangle
        //
        for (i = 1; i <= nb; i = i + 1) {
            //
            //           Update A(i:n,i)
            //
            Rgemv("No transpose", n - i + 1, i - 1, -one, a[(i - 1)], lda, w[(i - 1)], ldw, one, a[(i - 1) + (i - 1) * lda], 1);
            Rgemv("No transpose", n - i + 1, i - 1, -one, w[(i - 1)], ldw, a[(i - 1)], lda, one, a[(i - 1) + (i - 1) * lda], 1);
            if (i < n) {
                //
                //              Generate elementary reflector H(i) to annihilate
                //              A(i+2:n,i)
                //
                Rlarfg(n - i, a[((i + 1) - 1) + (i - 1) * lda], a[((min(i + 2) - 1) + (n)-1) * lda], 1, tau[i - 1]);
                e[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                a[((i + 1) - 1) + (i - 1) * lda] = one;
                //
                //              Compute W(i+1:n,i)
                //
                Rsymv("Lower", n - i, one, a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda, a[((i + 1) - 1) + (i - 1) * lda], 1, zero, w[((i + 1) - 1) + (i - 1) * ldw], 1);
                Rgemv("Transpose", n - i, i - 1, one, w[((i + 1) - 1)], ldw, a[((i + 1) - 1) + (i - 1) * lda], 1, zero, w[(i - 1) * ldw], 1);
                Rgemv("No transpose", n - i, i - 1, -one, a[((i + 1) - 1)], lda, w[(i - 1) * ldw], 1, one, w[((i + 1) - 1) + (i - 1) * ldw], 1);
                Rgemv("Transpose", n - i, i - 1, one, a[((i + 1) - 1)], lda, a[((i + 1) - 1) + (i - 1) * lda], 1, zero, w[(i - 1) * ldw], 1);
                Rgemv("No transpose", n - i, i - 1, -one, w[((i + 1) - 1)], ldw, w[(i - 1) * ldw], 1, one, w[((i + 1) - 1) + (i - 1) * ldw], 1);
                Rscal(n - i, tau[i - 1], w[((i + 1) - 1) + (i - 1) * ldw], 1);
                alpha = -half * tau[i - 1] * Rdot(n - i, w[((i + 1) - 1) + (i - 1) * ldw], 1, a[((i + 1) - 1) + (i - 1) * lda], 1);
                Raxpy(n - i, alpha, a[((i + 1) - 1) + (i - 1) * lda], 1, w[((i + 1) - 1) + (i - 1) * ldw], 1);
            }
            //
        }
    }
    //
    //     End of Rlatrd
    //
}
