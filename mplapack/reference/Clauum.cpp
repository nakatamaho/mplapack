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

void Clauum(const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER const &lda, INTEGER &info) {
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
        Mxerbla("Clauum", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Determine the block size for this environment.
    //
    INTEGER nb = iMlaenv[("Clauum" - 1) * ldiMlaenv];
    //
    INTEGER i = 0;
    INTEGER ib = 0;
    const COMPLEX cone = (1.0, 0.0);
    const REAL one = 1.0;
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code
        //
        Clauu2(uplo, n, a, lda, info);
    } else {
        //
        //        Use blocked code
        //
        if (upper) {
            //
            //           Compute the product U * U**H.
            //
            for (i = 1; i <= n; i = i + nb) {
                ib = min(nb, n - i + 1);
                Ctrmm("Right", "Upper", "Conjugate transpose", "Non-unit", i - 1, ib, cone, a[(i - 1) + (i - 1) * lda], lda, a[(i - 1) * lda], lda);
                Clauu2("Upper", ib, a[(i - 1) + (i - 1) * lda], lda, info);
                if (i + ib <= n) {
                    Cgemm("No transpose", "Conjugate transpose", i - 1, ib, n - i - ib + 1, cone, a[((i + ib) - 1) * lda], lda, a[(i - 1) + ((i + ib) - 1) * lda], lda, cone, a[(i - 1) * lda], lda);
                    Cherk("Upper", "No transpose", ib, n - i - ib + 1, one, a[(i - 1) + ((i + ib) - 1) * lda], lda, one, a[(i - 1) + (i - 1) * lda], lda);
                }
            }
        } else {
            //
            //           Compute the product L**H * L.
            //
            for (i = 1; i <= n; i = i + nb) {
                ib = min(nb, n - i + 1);
                Ctrmm("Left", "Lower", "Conjugate transpose", "Non-unit", ib, i - 1, cone, a[(i - 1) + (i - 1) * lda], lda, a[(i - 1)], lda);
                Clauu2("Lower", ib, a[(i - 1) + (i - 1) * lda], lda, info);
                if (i + ib <= n) {
                    Cgemm("Conjugate transpose", "No transpose", ib, i - 1, n - i - ib + 1, cone, a[((i + ib) - 1) + (i - 1) * lda], lda, a[((i + ib) - 1)], lda, cone, a[(i - 1)], lda);
                    Cherk("Lower", "Conjugate transpose", ib, n - i - ib + 1, one, a[((i + ib) - 1) + (i - 1) * lda], lda, one, a[(i - 1) + (i - 1) * lda], lda);
                }
            }
        }
    }
    //
    //     End of Clauum
    //
}
