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

void Rlauu2(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, INTEGER &info) {
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
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Rlauu2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    INTEGER i = 0;
    REAL aii = 0.0;
    const REAL one = 1.0;
    if (upper) {
        //
        //        Compute the product U * U**T.
        //
        for (i = 1; i <= n; i = i + 1) {
            aii = a[(i - 1) + (i - 1) * lda];
            if (i < n) {
                a[(i - 1) + (i - 1) * lda] = Rdot(n - i + 1, a[(i - 1) + (i - 1) * lda], lda, a[(i - 1) + (i - 1) * lda], lda);
                Rgemv("No transpose", i - 1, n - i, one, a[((i + 1) - 1) * lda], lda, a[(i - 1) + ((i + 1) - 1) * lda], lda, aii, a[(i - 1) * lda], 1);
            } else {
                Rscal(i, aii, a[(i - 1) * lda], 1);
            }
        }
        //
    } else {
        //
        //        Compute the product L**T * L.
        //
        for (i = 1; i <= n; i = i + 1) {
            aii = a[(i - 1) + (i - 1) * lda];
            if (i < n) {
                a[(i - 1) + (i - 1) * lda] = Rdot(n - i + 1, a[(i - 1) + (i - 1) * lda], 1, a[(i - 1) + (i - 1) * lda], 1);
                Rgemv("Transpose", n - i, i - 1, one, a[((i + 1) - 1)], lda, a[((i + 1) - 1) + (i - 1) * lda], 1, aii, a[(i - 1)], lda);
            } else {
                Rscal(i, aii, a[(i - 1)], lda);
            }
        }
    }
    //
    //     End of Rlauu2
    //
}
