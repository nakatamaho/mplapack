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

void Chegs2(INTEGER const &itype, const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER const &lda, COMPLEX *b, INTEGER const &ldb, INTEGER &info) {
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
        Mxerbla("Chegs2", -info);
        return;
    }
    //
    INTEGER k = 0;
    REAL akk = 0.0;
    REAL bkk = 0.0;
    const REAL one = 1.0;
    const REAL half = 0.5e+0;
    COMPLEX ct = 0.0;
    const COMPLEX cone = (1.0, 0.0);
    if (itype == 1) {
        if (upper) {
            //
            //           Compute inv(U**H)*A*inv(U)
            //
            for (k = 1; k <= n; k = k + 1) {
                //
                //              Update the upper triangle of A(k:n,k:n)
                //
                akk = a[(k - 1) + (k - 1) * lda];
                bkk = b[(k - 1) + (k - 1) * ldb];
                akk = akk / pow2(bkk);
                a[(k - 1) + (k - 1) * lda] = akk;
                if (k < n) {
                    CRscal(n - k, one / bkk, a[(k - 1) + ((k + 1) - 1) * lda], lda);
                    ct = -half * akk;
                    Clacgv(n - k, a[(k - 1) + ((k + 1) - 1) * lda], lda);
                    Clacgv(n - k, b[(k - 1) + ((k + 1) - 1) * ldb], ldb);
                    Caxpy(n - k, ct, b[(k - 1) + ((k + 1) - 1) * ldb], ldb, a[(k - 1) + ((k + 1) - 1) * lda], lda);
                    Cher2(uplo, n - k, -cone, a[(k - 1) + ((k + 1) - 1) * lda], lda, b[(k - 1) + ((k + 1) - 1) * ldb], ldb, a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda);
                    Caxpy(n - k, ct, b[(k - 1) + ((k + 1) - 1) * ldb], ldb, a[(k - 1) + ((k + 1) - 1) * lda], lda);
                    Clacgv(n - k, b[(k - 1) + ((k + 1) - 1) * ldb], ldb);
                    Ctrsv(uplo, "Conjugate transpose", "Non-unit", n - k, b[((k + 1) - 1) + ((k + 1) - 1) * ldb], ldb, a[(k - 1) + ((k + 1) - 1) * lda], lda);
                    Clacgv(n - k, a[(k - 1) + ((k + 1) - 1) * lda], lda);
                }
            }
        } else {
            //
            //           Compute inv(L)*A*inv(L**H)
            //
            for (k = 1; k <= n; k = k + 1) {
                //
                //              Update the lower triangle of A(k:n,k:n)
                //
                akk = a[(k - 1) + (k - 1) * lda];
                bkk = b[(k - 1) + (k - 1) * ldb];
                akk = akk / pow2(bkk);
                a[(k - 1) + (k - 1) * lda] = akk;
                if (k < n) {
                    CRscal(n - k, one / bkk, a[((k + 1) - 1) + (k - 1) * lda], 1);
                    ct = -half * akk;
                    Caxpy(n - k, ct, b[((k + 1) - 1) + (k - 1) * ldb], 1, a[((k + 1) - 1) + (k - 1) * lda], 1);
                    Cher2(uplo, n - k, -cone, a[((k + 1) - 1) + (k - 1) * lda], 1, b[((k + 1) - 1) + (k - 1) * ldb], 1, a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda);
                    Caxpy(n - k, ct, b[((k + 1) - 1) + (k - 1) * ldb], 1, a[((k + 1) - 1) + (k - 1) * lda], 1);
                    Ctrsv(uplo, "No transpose", "Non-unit", n - k, b[((k + 1) - 1) + ((k + 1) - 1) * ldb], ldb, a[((k + 1) - 1) + (k - 1) * lda], 1);
                }
            }
        }
    } else {
        if (upper) {
            //
            //           Compute U*A*U**H
            //
            for (k = 1; k <= n; k = k + 1) {
                //
                //              Update the upper triangle of A(1:k,1:k)
                //
                akk = a[(k - 1) + (k - 1) * lda];
                bkk = b[(k - 1) + (k - 1) * ldb];
                Ctrmv(uplo, "No transpose", "Non-unit", k - 1, b, ldb, a[(k - 1) * lda], 1);
                ct = half * akk;
                Caxpy(k - 1, ct, b[(k - 1) * ldb], 1, a[(k - 1) * lda], 1);
                Cher2(uplo, k - 1, cone, a[(k - 1) * lda], 1, b[(k - 1) * ldb], 1, a, lda);
                Caxpy(k - 1, ct, b[(k - 1) * ldb], 1, a[(k - 1) * lda], 1);
                CRscal(k - 1, bkk, a[(k - 1) * lda], 1);
                a[(k - 1) + (k - 1) * lda] = akk * pow2(bkk);
            }
        } else {
            //
            //           Compute L**H *A*L
            //
            for (k = 1; k <= n; k = k + 1) {
                //
                //              Update the lower triangle of A(1:k,1:k)
                //
                akk = a[(k - 1) + (k - 1) * lda];
                bkk = b[(k - 1) + (k - 1) * ldb];
                Clacgv(k - 1, a[(k - 1)], lda);
                Ctrmv(uplo, "Conjugate transpose", "Non-unit", k - 1, b, ldb, a[(k - 1)], lda);
                ct = half * akk;
                Clacgv(k - 1, b[(k - 1)], ldb);
                Caxpy(k - 1, ct, b[(k - 1)], ldb, a[(k - 1)], lda);
                Cher2(uplo, k - 1, cone, a[(k - 1)], lda, b[(k - 1)], ldb, a, lda);
                Caxpy(k - 1, ct, b[(k - 1)], ldb, a[(k - 1)], lda);
                Clacgv(k - 1, b[(k - 1)], ldb);
                CRscal(k - 1, bkk, a[(k - 1)], lda);
                Clacgv(k - 1, a[(k - 1)], lda);
                a[(k - 1) + (k - 1) * lda] = akk * pow2(bkk);
            }
        }
    }
    //
    //     End of Chegs2
    //
}
