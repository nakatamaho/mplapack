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

void Rpotrf2(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, INTEGER &info) {
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
    //     Test the input parameters
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
        Mxerbla("Rpotrf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     N=1 case
    //
    const REAL zero = 0.0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER iinfo = 0;
    const REAL one = 1.0;
    if (n == 1) {
        //
        //        Test for non-positive-definiteness
        //
        if (a[(1 - 1)] <= zero || disnan(a[(1 - 1)])) {
            info = 1;
            return;
        }
        //
        //        Factor
        //
        a[(1 - 1)] = sqrt(a[(1 - 1)]);
        //
        //     Use recursive code
        //
    } else {
        n1 = n / 2;
        n2 = n - n1;
        //
        //        Factor A11
        //
        Rpotrf2(uplo, n1, a[(1 - 1)], lda, iinfo);
        if (iinfo != 0) {
            info = iinfo;
            return;
        }
        //
        //        Compute the Cholesky factorization A = U**T*U
        //
        if (upper) {
            //
            //           Update and scale A12
            //
            Rtrsm("L", "U", "T", "N", n1, n2, one, a[(1 - 1)], lda, a[((n1 + 1) - 1) * lda], lda);
            //
            //           Update and factor A22
            //
            Rsyrk(uplo, "T", n2, n1, -one, a[((n1 + 1) - 1) * lda], lda, one, a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda);
            Rpotrf2(uplo, n2, a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda, iinfo);
            if (iinfo != 0) {
                info = iinfo + n1;
                return;
            }
            //
            //        Compute the Cholesky factorization A = L*L**T
            //
        } else {
            //
            //           Update and scale A21
            //
            Rtrsm("R", "L", "T", "N", n2, n1, one, a[(1 - 1)], lda, a[((n1 + 1) - 1)], lda);
            //
            //           Update and factor A22
            //
            Rsyrk(uplo, "N", n2, n1, -one, a[((n1 + 1) - 1)], lda, one, a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda);
            Rpotrf2(uplo, n2, a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda, iinfo);
            if (iinfo != 0) {
                info = iinfo + n1;
                return;
            }
        }
    }
    //
    //     End of Rpotrf2
    //
}
