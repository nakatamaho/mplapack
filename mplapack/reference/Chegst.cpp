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

void Chegst(INTEGER const itype, const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    if (itype < 1 || itype > 3) {
        info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Chegst", -info);
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
    INTEGER nb = iMlaenv(1, "Chegst", uplo, n, -1, -1, -1);
    //
    INTEGER k = 0;
    INTEGER kb = 0;
    const COMPLEX cone = (1.0, 0.0);
    const COMPLEX half = (0.5e+0, 0.0);
    const REAL one = 1.0;
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code
        //
        Chegs2(itype, uplo, n, a, lda, b, ldb, info);
    } else {
        //
        //        Use blocked code
        //
        if (itype == 1) {
            if (upper) {
                //
                //              Compute inv(U**H)*A*inv(U)
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the upper triangle of A(k:n,k:n)
                    //
                    Chegs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                    if (k + kb <= n) {
                        Ctrsm("Left", uplo, "Conjugate transpose", "Non-unit", kb, n - k - kb + 1, cone, &b[(k - 1) + (k - 1) * ldb], ldb, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                        Chemm("Left", uplo, kb, n - k - kb + 1, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + ((k + kb) - 1) * ldb], ldb, cone, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                        Cher2k(uplo, "Conjugate transpose", n - k - kb + 1, kb, -cone, &a[(k - 1) + ((k + kb) - 1) * lda], lda, &b[(k - 1) + ((k + kb) - 1) * ldb], ldb, one, &a[((k + kb) - 1) + ((k + kb) - 1) * lda], lda);
                        Chemm("Left", uplo, kb, n - k - kb + 1, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + ((k + kb) - 1) * ldb], ldb, cone, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                        Ctrsm("Right", uplo, "No transpose", "Non-unit", kb, n - k - kb + 1, cone, &b[((k + kb) - 1) + ((k + kb) - 1) * ldb], ldb, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                    }
                }
            } else {
                //
                //              Compute inv(L)*A*inv(L**H)
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the lower triangle of A(k:n,k:n)
                    //
                    Chegs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                    if (k + kb <= n) {
                        Ctrsm("Right", uplo, "Conjugate transpose", "Non-unit", n - k - kb + 1, kb, cone, &b[(k - 1) + (k - 1) * ldb], ldb, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                        Chemm("Right", uplo, n - k - kb + 1, kb, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[((k + kb) - 1) + (k - 1) * ldb], ldb, cone, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                        Cher2k(uplo, "No transpose", n - k - kb + 1, kb, -cone, &a[((k + kb) - 1) + (k - 1) * lda], lda, &b[((k + kb) - 1) + (k - 1) * ldb], ldb, one, &a[((k + kb) - 1) + ((k + kb) - 1) * lda], lda);
                        Chemm("Right", uplo, n - k - kb + 1, kb, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[((k + kb) - 1) + (k - 1) * ldb], ldb, cone, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                        Ctrsm("Left", uplo, "No transpose", "Non-unit", n - k - kb + 1, kb, cone, &b[((k + kb) - 1) + ((k + kb) - 1) * ldb], ldb, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                    }
                }
            }
        } else {
            if (upper) {
                //
                //              Compute U*A*U**H
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
                    //
                    Ctrmm("Left", uplo, "No transpose", "Non-unit", k - 1, kb, cone, b, ldb, &a[(k - 1) * lda], lda);
                    Chemm("Right", uplo, k - 1, kb, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) * ldb], ldb, cone, &a[(k - 1) * lda], lda);
                    Cher2k(uplo, "No transpose", k - 1, kb, cone, &a[(k - 1) * lda], lda, &b[(k - 1) * ldb], ldb, one, a, lda);
                    Chemm("Right", uplo, k - 1, kb, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) * ldb], ldb, cone, &a[(k - 1) * lda], lda);
                    Ctrmm("Right", uplo, "Conjugate transpose", "Non-unit", k - 1, kb, cone, &b[(k - 1) + (k - 1) * ldb], ldb, &a[(k - 1) * lda], lda);
                    Chegs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                }
            } else {
                //
                //              Compute L**H*A*L
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
                    //
                    Ctrmm("Right", uplo, "No transpose", "Non-unit", kb, k - 1, cone, b, ldb, &a[(k - 1)], lda);
                    Chemm("Left", uplo, kb, k - 1, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1)], ldb, cone, &a[(k - 1)], lda);
                    Cher2k(uplo, "Conjugate transpose", k - 1, kb, cone, &a[(k - 1)], lda, &b[(k - 1)], ldb, one, a, lda);
                    Chemm("Left", uplo, kb, k - 1, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1)], ldb, cone, &a[(k - 1)], lda);
                    Ctrmm("Left", uplo, "Conjugate transpose", "Non-unit", kb, k - 1, cone, &b[(k - 1) + (k - 1) * ldb], ldb, &a[(k - 1)], lda);
                    Chegs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                }
            }
        }
    }
    //
    //     End of Chegst
    //
}
