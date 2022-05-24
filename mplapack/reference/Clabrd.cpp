/*
 * Copyright (c) 2008-2022
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

void Clabrd(INTEGER const m, INTEGER const n, INTEGER const nb, COMPLEX *a, INTEGER const lda, REAL *d, REAL *e, COMPLEX *tauq, COMPLEX *taup, COMPLEX *x, INTEGER const ldx, COMPLEX *y, INTEGER const ldy) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    INTEGER i = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    COMPLEX alpha = 0.0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    if (m >= n) {
        //
        //        Reduce to upper bidiagonal form
        //
        for (i = 1; i <= nb; i = i + 1) {
            //
            //           Update A(i:m,i)
            //
            Clacgv(i - 1, &y[(i - 1)], ldy);
            Cgemv("No transpose", m - i + 1, i - 1, -one, &a[(i - 1)], lda, &y[(i - 1)], ldy, one, &a[(i - 1) + (i - 1) * lda], 1);
            Clacgv(i - 1, &y[(i - 1)], ldy);
            Cgemv("No transpose", m - i + 1, i - 1, -one, &x[(i - 1)], ldx, &a[(i - 1) * lda], 1, one, &a[(i - 1) + (i - 1) * lda], 1);
            //
            //           Generate reflection Q(i) to annihilate A(i+1:m,i)
            //
            alpha = a[(i - 1) + (i - 1) * lda];
            Clarfg(m - i + 1, alpha, &a[(min(i + 1, m) - 1) + (i - 1) * lda], 1, tauq[i - 1]);
            d[i - 1] = alpha.real();
            if (i < n) {
                a[(i - 1) + (i - 1) * lda] = one;
                //
                //              Compute Y(i+1:n,i)
                //
                Cgemv("Conjugate transpose", m - i + 1, n - i, one, &a[(i - 1) + ((i + 1) - 1) * lda], lda, &a[(i - 1) + (i - 1) * lda], 1, zero, &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                Cgemv("Conjugate transpose", m - i + 1, i - 1, one, &a[(i - 1)], lda, &a[(i - 1) + (i - 1) * lda], 1, zero, &y[(i - 1) * ldy], 1);
                Cgemv("No transpose", n - i, i - 1, -one, &y[((i + 1) - 1)], ldy, &y[(i - 1) * ldy], 1, one, &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                Cgemv("Conjugate transpose", m - i + 1, i - 1, one, &x[(i - 1)], ldx, &a[(i - 1) + (i - 1) * lda], 1, zero, &y[(i - 1) * ldy], 1);
                Cgemv("Conjugate transpose", i - 1, n - i, -one, &a[((i + 1) - 1) * lda], lda, &y[(i - 1) * ldy], 1, one, &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                Cscal(n - i, tauq[i - 1], &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                //
                //              Update A(i,i+1:n)
                //
                Clacgv(n - i, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
                Clacgv(i, &a[(i - 1)], lda);
                Cgemv("No transpose", n - i, i, -one, &y[((i + 1) - 1)], ldy, &a[(i - 1)], lda, one, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
                Clacgv(i, &a[(i - 1)], lda);
                Clacgv(i - 1, &x[(i - 1)], ldx);
                Cgemv("Conjugate transpose", i - 1, n - i, -one, &a[((i + 1) - 1) * lda], lda, &x[(i - 1)], ldx, one, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
                Clacgv(i - 1, &x[(i - 1)], ldx);
                //
                //              Generate reflection P(i) to annihilate A(i,i+2:n)
                //
                alpha = a[(i - 1) + ((i + 1) - 1) * lda];
                Clarfg(n - i, alpha, &a[(i - 1) + (min(i + 2, n) - 1) * lda], lda, taup[i - 1]);
                e[i - 1] = alpha.real();
                a[(i - 1) + ((i + 1) - 1) * lda] = one;
                //
                //              Compute X(i+1:m,i)
                //
                Cgemv("No transpose", m - i, n - i, one, &a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda, &a[(i - 1) + ((i + 1) - 1) * lda], lda, zero, &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Cgemv("Conjugate transpose", n - i, i, one, &y[((i + 1) - 1)], ldy, &a[(i - 1) + ((i + 1) - 1) * lda], lda, zero, &x[(i - 1) * ldx], 1);
                Cgemv("No transpose", m - i, i, -one, &a[((i + 1) - 1)], lda, &x[(i - 1) * ldx], 1, one, &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Cgemv("No transpose", i - 1, n - i, one, &a[((i + 1) - 1) * lda], lda, &a[(i - 1) + ((i + 1) - 1) * lda], lda, zero, &x[(i - 1) * ldx], 1);
                Cgemv("No transpose", m - i, i - 1, -one, &x[((i + 1) - 1)], ldx, &x[(i - 1) * ldx], 1, one, &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Cscal(m - i, taup[i - 1], &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Clacgv(n - i, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
            }
        }
    } else {
        //
        //        Reduce to lower bidiagonal form
        //
        for (i = 1; i <= nb; i = i + 1) {
            //
            //           Update A(i,i:n)
            //
            Clacgv(n - i + 1, &a[(i - 1) + (i - 1) * lda], lda);
            Clacgv(i - 1, &a[(i - 1)], lda);
            Cgemv("No transpose", n - i + 1, i - 1, -one, &y[(i - 1)], ldy, &a[(i - 1)], lda, one, &a[(i - 1) + (i - 1) * lda], lda);
            Clacgv(i - 1, &a[(i - 1)], lda);
            Clacgv(i - 1, &x[(i - 1)], ldx);
            Cgemv("Conjugate transpose", i - 1, n - i + 1, -one, &a[(i - 1) * lda], lda, &x[(i - 1)], ldx, one, &a[(i - 1) + (i - 1) * lda], lda);
            Clacgv(i - 1, &x[(i - 1)], ldx);
            //
            //           Generate reflection P(i) to annihilate A(i,i+1:n)
            //
            alpha = a[(i - 1) + (i - 1) * lda];
            Clarfg(n - i + 1, alpha, &a[(i - 1) + (min(i + 1, n) - 1) * lda], lda, taup[i - 1]);
            d[i - 1] = alpha.real();
            if (i < m) {
                a[(i - 1) + (i - 1) * lda] = one;
                //
                //              Compute X(i+1:m,i)
                //
                Cgemv("No transpose", m - i, n - i + 1, one, &a[((i + 1) - 1) + (i - 1) * lda], lda, &a[(i - 1) + (i - 1) * lda], lda, zero, &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Cgemv("Conjugate transpose", n - i + 1, i - 1, one, &y[(i - 1)], ldy, &a[(i - 1) + (i - 1) * lda], lda, zero, &x[(i - 1) * ldx], 1);
                Cgemv("No transpose", m - i, i - 1, -one, &a[((i + 1) - 1)], lda, &x[(i - 1) * ldx], 1, one, &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Cgemv("No transpose", i - 1, n - i + 1, one, &a[(i - 1) * lda], lda, &a[(i - 1) + (i - 1) * lda], lda, zero, &x[(i - 1) * ldx], 1);
                Cgemv("No transpose", m - i, i - 1, -one, &x[((i + 1) - 1)], ldx, &x[(i - 1) * ldx], 1, one, &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Cscal(m - i, taup[i - 1], &x[((i + 1) - 1) + (i - 1) * ldx], 1);
                Clacgv(n - i + 1, &a[(i - 1) + (i - 1) * lda], lda);
                //
                //              Update A(i+1:m,i)
                //
                Clacgv(i - 1, &y[(i - 1)], ldy);
                Cgemv("No transpose", m - i, i - 1, -one, &a[((i + 1) - 1)], lda, &y[(i - 1)], ldy, one, &a[((i + 1) - 1) + (i - 1) * lda], 1);
                Clacgv(i - 1, &y[(i - 1)], ldy);
                Cgemv("No transpose", m - i, i, -one, &x[((i + 1) - 1)], ldx, &a[(i - 1) * lda], 1, one, &a[((i + 1) - 1) + (i - 1) * lda], 1);
                //
                //              Generate reflection Q(i) to annihilate A(i+2:m,i)
                //
                alpha = a[((i + 1) - 1) + (i - 1) * lda];
                Clarfg(m - i, alpha, &a[(min(i + 2, m) - 1) + (i - 1) * lda], 1, tauq[i - 1]);
                e[i - 1] = alpha.real();
                a[((i + 1) - 1) + (i - 1) * lda] = one;
                //
                //              Compute Y(i+1:n,i)
                //
                Cgemv("Conjugate transpose", m - i, n - i, one, &a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda, &a[((i + 1) - 1) + (i - 1) * lda], 1, zero, &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                Cgemv("Conjugate transpose", m - i, i - 1, one, &a[((i + 1) - 1)], lda, &a[((i + 1) - 1) + (i - 1) * lda], 1, zero, &y[(i - 1) * ldy], 1);
                Cgemv("No transpose", n - i, i - 1, -one, &y[((i + 1) - 1)], ldy, &y[(i - 1) * ldy], 1, one, &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                Cgemv("Conjugate transpose", m - i, i, one, &x[((i + 1) - 1)], ldx, &a[((i + 1) - 1) + (i - 1) * lda], 1, zero, &y[(i - 1) * ldy], 1);
                Cgemv("Conjugate transpose", i, n - i, -one, &a[((i + 1) - 1) * lda], lda, &y[(i - 1) * ldy], 1, one, &y[((i + 1) - 1) + (i - 1) * ldy], 1);
                Cscal(n - i, tauq[i - 1], &y[((i + 1) - 1) + (i - 1) * ldy], 1);
            } else {
                Clacgv(n - i + 1, &a[(i - 1) + (i - 1) * lda], lda);
            }
        }
    }
    //
    //     End of Clabrd
    //
}
