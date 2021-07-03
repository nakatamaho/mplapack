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

void Rsygst(INTEGER const itype, const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER &info) {
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
        Mxerbla("Rsygst", -info);
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
    INTEGER nb = iMlaenv(1, "Rsygst", uplo, n, -1, -1, -1);
    //
    INTEGER k = 0;
    INTEGER kb = 0;
    const REAL one = 1.0;
    const REAL half = 0.5e0;
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code
        //
        Rsygs2(itype, uplo, n, a, lda, b, ldb, info);
    } else {
        //
        //        Use blocked code
        //
        if (itype == 1) {
            if (upper) {
                //
                //              Compute inv(U**T)*A*inv(U)
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the upper triangle of A(k:n,k:n)
                    //
                    Rsygs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                    if (k + kb <= n) {
                        Rtrsm("Left", uplo, "Transpose", "Non-unit", kb, n - k - kb + 1, one, &b[(k - 1) + (k - 1) * ldb], ldb, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                        Rsymm("Left", uplo, kb, n - k - kb + 1, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + ((k + kb) - 1) * ldb], ldb, one, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                        Rsyr2k(uplo, "Transpose", n - k - kb + 1, kb, -one, &a[(k - 1) + ((k + kb) - 1) * lda], lda, &b[(k - 1) + ((k + kb) - 1) * ldb], ldb, one, &a[((k + kb) - 1) + ((k + kb) - 1) * lda], lda);
                        Rsymm("Left", uplo, kb, n - k - kb + 1, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + ((k + kb) - 1) * ldb], ldb, one, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                        Rtrsm("Right", uplo, "No transpose", "Non-unit", kb, n - k - kb + 1, one, &b[((k + kb) - 1) + ((k + kb) - 1) * ldb], ldb, &a[(k - 1) + ((k + kb) - 1) * lda], lda);
                    }
                }
            } else {
                //
                //              Compute inv(L)*A*inv(L**T)
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the lower triangle of A(k:n,k:n)
                    //
                    Rsygs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                    if (k + kb <= n) {
                        Rtrsm("Right", uplo, "Transpose", "Non-unit", n - k - kb + 1, kb, one, &b[(k - 1) + (k - 1) * ldb], ldb, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                        Rsymm("Right", uplo, n - k - kb + 1, kb, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[((k + kb) - 1) + (k - 1) * ldb], ldb, one, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                        Rsyr2k(uplo, "No transpose", n - k - kb + 1, kb, -one, &a[((k + kb) - 1) + (k - 1) * lda], lda, &b[((k + kb) - 1) + (k - 1) * ldb], ldb, one, &a[((k + kb) - 1) + ((k + kb) - 1) * lda], lda);
                        Rsymm("Right", uplo, n - k - kb + 1, kb, -half, &a[(k - 1) + (k - 1) * lda], lda, &b[((k + kb) - 1) + (k - 1) * ldb], ldb, one, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                        Rtrsm("Left", uplo, "No transpose", "Non-unit", n - k - kb + 1, kb, one, &b[((k + kb) - 1) + ((k + kb) - 1) * ldb], ldb, &a[((k + kb) - 1) + (k - 1) * lda], lda);
                    }
                }
            }
        } else {
            if (upper) {
                //
                //              Compute U*A*U**T
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
                    //
                    Rtrmm("Left", uplo, "No transpose", "Non-unit", k - 1, kb, one, b, ldb, &a[(k - 1) * lda], lda);
                    Rsymm("Right", uplo, k - 1, kb, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) * ldb], ldb, one, &a[(k - 1) * lda], lda);
                    Rsyr2k(uplo, "No transpose", k - 1, kb, one, &a[(k - 1) * lda], lda, &b[(k - 1) * ldb], ldb, one, a, lda);
                    Rsymm("Right", uplo, k - 1, kb, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) * ldb], ldb, one, &a[(k - 1) * lda], lda);
                    Rtrmm("Right", uplo, "Transpose", "Non-unit", k - 1, kb, one, &b[(k - 1) + (k - 1) * ldb], ldb, &a[(k - 1) * lda], lda);
                    Rsygs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                }
            } else {
                //
                //              Compute L**T*A*L
                //
                for (k = 1; k <= n; k = k + nb) {
                    kb = min(n - k + 1, nb);
                    //
                    //                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
                    //
                    Rtrmm("Right", uplo, "No transpose", "Non-unit", kb, k - 1, one, b, ldb, &a[(k - 1)], lda);
                    Rsymm("Left", uplo, kb, k - 1, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1)], ldb, one, &a[(k - 1)], lda);
                    Rsyr2k(uplo, "Transpose", k - 1, kb, one, &a[(k - 1)], lda, &b[(k - 1)], ldb, one, a, lda);
                    Rsymm("Left", uplo, kb, k - 1, half, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1)], ldb, one, &a[(k - 1)], lda);
                    Rtrmm("Left", uplo, "Transpose", "Non-unit", kb, k - 1, one, &b[(k - 1) + (k - 1) * ldb], ldb, &a[(k - 1)], lda);
                    Rsygs2(itype, uplo, kb, &a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, info);
                }
            }
        }
    }
    //
    //     End of Rsygst
    //
}
