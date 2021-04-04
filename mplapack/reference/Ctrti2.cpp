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

void Ctrti2(const char *uplo, const char *diag, INTEGER const &n, COMPLEX *a, INTEGER const &lda, INTEGER &info) {
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
    bool upper = Mlsame(uplo, "U");
    bool nounit = Mlsame(diag, "N");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (!nounit && !Mlsame(diag, "U")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Ctrti2", -info);
        return;
    }
    //
    INTEGER j = 0;
    const COMPLEX one = (1.0, 0.0);
    COMPLEX ajj = 0.0;
    if (upper) {
        //
        //        Compute inverse of upper triangular matrix.
        //
        for (j = 1; j <= n; j = j + 1) {
            if (nounit) {
                a[(j - 1) + (j - 1) * lda] = one / a[(j - 1) + (j - 1) * lda];
                ajj = -a[(j - 1) + (j - 1) * lda];
            } else {
                ajj = -one;
            }
            //
            //           Compute elements 1:j-1 of j-th column.
            //
            Ctrmv("Upper", "No transpose", diag, j - 1, a, lda, a[(j - 1) * lda], 1);
            Cscal(j - 1, ajj, a[(j - 1) * lda], 1);
        }
    } else {
        //
        //        Compute inverse of lower triangular matrix.
        //
        for (j = n; j >= 1; j = j - 1) {
            if (nounit) {
                a[(j - 1) + (j - 1) * lda] = one / a[(j - 1) + (j - 1) * lda];
                ajj = -a[(j - 1) + (j - 1) * lda];
            } else {
                ajj = -one;
            }
            if (j < n) {
                //
                //              Compute elements j+1:n of j-th column.
                //
                Ctrmv("Lower", "No transpose", diag, n - j, a[((j + 1) - 1) + ((j + 1) - 1) * lda], lda, a[((j + 1) - 1) + (j - 1) * lda], 1);
                Cscal(n - j, ajj, a[((j + 1) - 1) + (j - 1) * lda], 1);
            }
        }
    }
    //
    //     End of Ctrti2
    //
}
