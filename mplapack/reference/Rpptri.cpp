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

void Rpptri(const char *uplo, INTEGER const &n, REAL *ap, INTEGER &info) {
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
    }
    if (info != 0) {
        Mxerbla("Rpptri", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Invert the triangular Cholesky factor U or L.
    //
    Rtptri(uplo, "Non-unit", n, ap, info);
    if (info > 0) {
        return;
    }
    //
    INTEGER jj = 0;
    INTEGER j = 0;
    INTEGER jc = 0;
    const REAL one = 1.0;
    REAL ajj = 0.0;
    INTEGER jjn = 0;
    if (upper) {
        //
        //        Compute the product inv(U) * inv(U)**T.
        //
        jj = 0;
        for (j = 1; j <= n; j = j + 1) {
            jc = jj + 1;
            jj += j;
            if (j > 1) {
                Rspr("Upper", j - 1, one, ap[jc - 1], 1, ap);
            }
            ajj = ap[jj - 1];
            Rscal(j, ajj, ap[jc - 1], 1);
        }
        //
    } else {
        //
        //        Compute the product inv(L)**T * inv(L).
        //
        jj = 1;
        for (j = 1; j <= n; j = j + 1) {
            jjn = jj + n - j + 1;
            ap[jj - 1] = Rdot(n - j + 1, ap[jj - 1], 1, ap[jj - 1], 1);
            if (j < n) {
                Rtpmv("Lower", "Transpose", "Non-unit", n - j, ap[jjn - 1], ap[(jj + 1) - 1], 1);
            }
            jj = jjn;
        }
    }
    //
    //     End of Rpptri
    //
}
