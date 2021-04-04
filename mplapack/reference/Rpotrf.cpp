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

void Rpotrf(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, INTEGER &info) {
    bool upper = false;
    INTEGER nb = 0;
    INTEGER j = 0;
    INTEGER jb = 0;
    const REAL one = 1.0;
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
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Rpotrf", -info);
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
    nb = iMlaenv[("Rpotrf" - 1) * ldiMlaenv];
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code.
        //
        Rpotrf2(uplo, n, a, lda, info);
    } else {
        //
        //        Use blocked code.
        //
        if (upper) {
            //
            //           Compute the Cholesky factorization A = U**T*U.
            //
            for (j = 1; j <= n; j = j + nb) {
                //
                //              Update and factorize the current diagonal block and test
                //              for non-positive-definiteness.
                //
                jb = min(nb, n - j + 1);
                Rsyrk("Upper", "Transpose", jb, j - 1, -one, a[(j - 1) * lda], lda, one, a[(j - 1) + (j - 1) * lda], lda);
                Rpotrf2("Upper", jb, a[(j - 1) + (j - 1) * lda], lda, info);
                if (info != 0) {
                    goto statement_30;
                }
                if (j + jb <= n) {
                    //
                    //                 Compute the current block row.
                    //
                    Rgemm("Transpose", "No transpose", jb, n - j - jb + 1, j - 1, -one, a[(j - 1) * lda], lda, a[((j + jb) - 1) * lda], lda, one, a[(j - 1) + ((j + jb) - 1) * lda], lda);
                    Rtrsm("Left", "Upper", "Transpose", "Non-unit", jb, n - j - jb + 1, one, a[(j - 1) + (j - 1) * lda], lda, a[(j - 1) + ((j + jb) - 1) * lda], lda);
                }
            }
            //
        } else {
            //
            //           Compute the Cholesky factorization A = L*L**T.
            //
            for (j = 1; j <= n; j = j + nb) {
                //
                //              Update and factorize the current diagonal block and test
                //              for non-positive-definiteness.
                //
                jb = min(nb, n - j + 1);
                Rsyrk("Lower", "No transpose", jb, j - 1, -one, a[(j - 1)], lda, one, a[(j - 1) + (j - 1) * lda], lda);
                Rpotrf2("Lower", jb, a[(j - 1) + (j - 1) * lda], lda, info);
                if (info != 0) {
                    goto statement_30;
                }
                if (j + jb <= n) {
                    //
                    //                 Compute the current block column.
                    //
                    Rgemm("No transpose", "Transpose", n - j - jb + 1, jb, j - 1, -one, a[((j + jb) - 1)], lda, a[(j - 1)], lda, one, a[((j + jb) - 1) + (j - 1) * lda], lda);
                    Rtrsm("Right", "Lower", "Transpose", "Non-unit", n - j - jb + 1, jb, one, a[(j - 1) + (j - 1) * lda], lda, a[((j + jb) - 1) + (j - 1) * lda], lda);
                }
            }
        }
    }
    goto statement_40;
//
statement_30:
    info += j - 1;
//
statement_40:;
    //
    //     End of Rpotrf
    //
}
