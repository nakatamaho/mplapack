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

void Csytrs_rook(const char *uplo, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, COMPLEX *b, INTEGER const ldb, INTEGER &info) {
    bool upper = false;
    INTEGER k = 0;
    INTEGER kp = 0;
    const COMPLEX cone = (1.0, 0.0);
    COMPLEX akm1k = 0.0;
    COMPLEX akm1 = 0.0;
    COMPLEX ak = 0.0;
    COMPLEX denom = 0.0;
    INTEGER j = 0;
    COMPLEX bkm1 = 0.0;
    COMPLEX bk = 0.0;
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
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Csytrs_rook", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    if (upper) {
        //
        //        Solve A*X = B, where A = U*D*U**T.
        //
        //        First solve U*D*X = B, overwriting B with X.
        //
        //        K is the main loop index, decreasing from N to 1 in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = n;
    statement_10:
        //
        //        If K < 1, exit from loop.
        //
        if (k < 1) {
            goto statement_30;
        }
        //
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Interchange rows K and IPIV(K).
            //
            kp = ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            //           Multiply by inv(U(K)), where U(K) is the transformation
            //           stored in column K of A.
            //
            Cgeru(k - 1, nrhs, -cone, &a[(k - 1) * lda], 1, &b[(k - 1)], ldb, &b[(1 - 1)], ldb);
            //
            //           Multiply by the inverse of the diagonal block.
            //
            Cscal(nrhs, cone / a[(k - 1) + (k - 1) * lda], &b[(k - 1)], ldb);
            k = k - 1;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)
            //
            kp = -ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            kp = -ipiv[(k - 1) - 1];
            if (kp != k - 1) {
                Cswap(nrhs, &b[((k - 1) - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            //           Multiply by inv(U(K)), where U(K) is the transformation
            //           stored in columns K-1 and K of A.
            //
            if (k > 2) {
                Cgeru(k - 2, nrhs, -cone, &a[(k - 1) * lda], 1, &b[(k - 1)], ldb, &b[(1 - 1)], ldb);
                Cgeru(k - 2, nrhs, -cone, &a[((k - 1) - 1) * lda], 1, &b[((k - 1) - 1)], ldb, &b[(1 - 1)], ldb);
            }
            //
            //           Multiply by the inverse of the diagonal block.
            //
            akm1k = a[((k - 1) - 1) + (k - 1) * lda];
            akm1 = a[((k - 1) - 1) + ((k - 1) - 1) * lda] / akm1k;
            ak = a[(k - 1) + (k - 1) * lda] / akm1k;
            denom = akm1 * ak - cone;
            for (j = 1; j <= nrhs; j = j + 1) {
                bkm1 = b[((k - 1) - 1) + (j - 1) * ldb] / akm1k;
                bk = b[(k - 1) + (j - 1) * ldb] / akm1k;
                b[((k - 1) - 1) + (j - 1) * ldb] = (ak * bkm1 - bk) / denom;
                b[(k - 1) + (j - 1) * ldb] = (akm1 * bk - bkm1) / denom;
            }
            k = k - 2;
        }
        //
        goto statement_10;
    statement_30:
        //
        //        Next solve U**T *X = B, overwriting B with X.
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = 1;
    statement_40:
        //
        //        If K > N, exit from loop.
        //
        if (k > n) {
            goto statement_50;
        }
        //
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Multiply by inv(U**T(K)), where U(K) is the transformation
            //           stored in column K of A.
            //
            if (k > 1) {
                Cgemv("Transpose", k - 1, nrhs, -cone, b, ldb, &a[(k - 1) * lda], 1, cone, &b[(k - 1)], ldb);
            }
            //
            //           Interchange rows K and IPIV(K).
            //
            kp = ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            k++;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
            //           stored in columns K and K+1 of A.
            //
            if (k > 1) {
                Cgemv("Transpose", k - 1, nrhs, -cone, b, ldb, &a[(k - 1) * lda], 1, cone, &b[(k - 1)], ldb);
                Cgemv("Transpose", k - 1, nrhs, -cone, b, ldb, &a[((k + 1) - 1) * lda], 1, cone, &b[((k + 1) - 1)], ldb);
            }
            //
            //           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1).
            //
            kp = -ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            kp = -ipiv[(k + 1) - 1];
            if (kp != k + 1) {
                Cswap(nrhs, &b[((k + 1) - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            k += 2;
        }
        //
        goto statement_40;
    statement_50:;
        //
    } else {
        //
        //        Solve A*X = B, where A = L*D*L**T.
        //
        //        First solve L*D*X = B, overwriting B with X.
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = 1;
    statement_60:
        //
        //        If K > N, exit from loop.
        //
        if (k > n) {
            goto statement_80;
        }
        //
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Interchange rows K and IPIV(K).
            //
            kp = ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            //           Multiply by inv(L(K)), where L(K) is the transformation
            //           stored in column K of A.
            //
            if (k < n) {
                Cgeru(n - k, nrhs, -cone, &a[((k + 1) - 1) + (k - 1) * lda], 1, &b[(k - 1)], ldb, &b[((k + 1) - 1)], ldb);
            }
            //
            //           Multiply by the inverse of the diagonal block.
            //
            Cscal(nrhs, cone / a[(k - 1) + (k - 1) * lda], &b[(k - 1)], ldb);
            k++;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1)
            //
            kp = -ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            kp = -ipiv[(k + 1) - 1];
            if (kp != k + 1) {
                Cswap(nrhs, &b[((k + 1) - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            //           Multiply by inv(L(K)), where L(K) is the transformation
            //           stored in columns K and K+1 of A.
            //
            if (k < n - 1) {
                Cgeru(n - k - 1, nrhs, -cone, &a[((k + 2) - 1) + (k - 1) * lda], 1, &b[(k - 1)], ldb, &b[((k + 2) - 1)], ldb);
                Cgeru(n - k - 1, nrhs, -cone, &a[((k + 2) - 1) + ((k + 1) - 1) * lda], 1, &b[((k + 1) - 1)], ldb, &b[((k + 2) - 1)], ldb);
            }
            //
            //           Multiply by the inverse of the diagonal block.
            //
            akm1k = a[((k + 1) - 1) + (k - 1) * lda];
            akm1 = a[(k - 1) + (k - 1) * lda] / akm1k;
            ak = a[((k + 1) - 1) + ((k + 1) - 1) * lda] / akm1k;
            denom = akm1 * ak - cone;
            for (j = 1; j <= nrhs; j = j + 1) {
                bkm1 = b[(k - 1) + (j - 1) * ldb] / akm1k;
                bk = b[((k + 1) - 1) + (j - 1) * ldb] / akm1k;
                b[(k - 1) + (j - 1) * ldb] = (ak * bkm1 - bk) / denom;
                b[((k + 1) - 1) + (j - 1) * ldb] = (akm1 * bk - bkm1) / denom;
            }
            k += 2;
        }
        //
        goto statement_60;
    statement_80:
        //
        //        Next solve L**T *X = B, overwriting B with X.
        //
        //        K is the main loop index, decreasing from N to 1 in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = n;
    statement_90:
        //
        //        If K < 1, exit from loop.
        //
        if (k < 1) {
            goto statement_100;
        }
        //
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Multiply by inv(L**T(K)), where L(K) is the transformation
            //           stored in column K of A.
            //
            if (k < n) {
                Cgemv("Transpose", n - k, nrhs, -cone, &b[((k + 1) - 1)], ldb, &a[((k + 1) - 1) + (k - 1) * lda], 1, cone, &b[(k - 1)], ldb);
            }
            //
            //           Interchange rows K and IPIV(K).
            //
            kp = ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            k = k - 1;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
            //           stored in columns K-1 and K of A.
            //
            if (k < n) {
                Cgemv("Transpose", n - k, nrhs, -cone, &b[((k + 1) - 1)], ldb, &a[((k + 1) - 1) + (k - 1) * lda], 1, cone, &b[(k - 1)], ldb);
                Cgemv("Transpose", n - k, nrhs, -cone, &b[((k + 1) - 1)], ldb, &a[((k + 1) - 1) + ((k - 1) - 1) * lda], 1, cone, &b[((k - 1) - 1)], ldb);
            }
            //
            //           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)
            //
            kp = -ipiv[k - 1];
            if (kp != k) {
                Cswap(nrhs, &b[(k - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            kp = -ipiv[(k - 1) - 1];
            if (kp != k - 1) {
                Cswap(nrhs, &b[((k - 1) - 1)], ldb, &b[(kp - 1)], ldb);
            }
            //
            k = k - 2;
        }
        //
        goto statement_90;
    statement_100:;
    }
    //
    //     End of Csytrs_rook
    //
}
