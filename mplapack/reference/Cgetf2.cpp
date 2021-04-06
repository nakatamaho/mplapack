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

void Cgetf2(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info) {
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
    //     Test the input parameters.
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
        Mxerbla("Cgetf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Compute machine safe minimum
    //
    REAL sfmin = Rlamch("S");
    //
    INTEGER j = 0;
    INTEGER jp = 0;
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    INTEGER i = 0;
    for (j = 1; j <= min(m, n); j = j + 1) {
        //
        //        Find pivot and test for singularity.
        //
        jp = j - 1 + iCamax(m - j + 1, &a[(j - 1) + (j - 1) * lda], 1);
        ipiv[j - 1] = jp;
        if (a[(jp - 1) + (j - 1) * lda] != zero) {
            //
            //           Apply the INTEGERerchange to columns 1:N.
            //
            if (jp != j) {
                Cswap(n, &a[(j - 1)], lda, &a[(jp - 1)], lda);
            }
            //
            //           Compute elements J+1:M of J-th column.
            //
            if (j < m) {
                if (abs(a[(j - 1) + (j - 1) * lda]) >= sfmin) {
                    Cscal(m - j, one / a[(j - 1) + (j - 1) * lda], &a[((j + 1) - 1) + (j - 1) * lda], 1);
                } else {
                    for (i = 1; i <= m - j; i = i + 1) {
                        a[((j + i) - 1) + (j - 1) * lda] = a[((j + i) - 1) + (j - 1) * lda] / a[(j - 1) + (j - 1) * lda];
                    }
                }
            }
            //
        } else if (info == 0) {
            //
            info = j;
        }
        //
        if (j < min(m, n)) {
            //
            //           Update trailing submatrix.
            //
            Cgeru(m - j, n - j, -one, &a[((j + 1) - 1) + (j - 1) * lda], 1, &a[(j - 1) + ((j + 1) - 1) * lda], lda, &a[((j + 1) - 1) + ((j + 1) - 1) * lda], lda);
        }
    }
    //
    //     End of Cgetf2
    //
}
