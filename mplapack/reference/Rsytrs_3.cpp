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

void Rsytrs_3(const char *uplo, INTEGER const &n, INTEGER const &nrhs, REAL *a, INTEGER const &lda, REAL *e, INTEGER *ipiv, REAL *b, INTEGER const &ldb, INTEGER &info) {
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
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rsytrs_3", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    INTEGER k = 0;
    INTEGER kp = 0;
    const REAL one = 1.0;
    INTEGER i = 0;
    REAL akm1k = 0.0;
    REAL akm1 = 0.0;
    REAL ak = 0.0;
    REAL denom = 0.0;
    INTEGER j = 0;
    REAL bkm1 = 0.0;
    REAL bk = 0.0;
    if (upper) {
        //
        //        Begin Upper
        //
        //        Solve A*X = B, where A = U*D*U**T.
        //
        //        P**T * B
        //
        //        Interchange rows K and IPIV(K) of matrix B in the same order
        //        that the formation order of IPIV(I) vector for Upper case.
        //
        //        (We can do the simple loop over IPIV with decrement -1,
        //        since the ABS value of IPIV( I ) represents the row index
        //        of the INTEGERerchange with row i in both 1x1 and 2x2 pivot cases)
        //
        for (k = n; k >= 1; k = k - 1) {
            kp = abs(ipiv[k - 1]);
            if (kp != k) {
                Rswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
            }
        }
        //
        //        Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
        //
        Rtrsm("L", "U", "N", "U", n, nrhs, one, a, lda, b, ldb);
        //
        //        Compute D \ B -> B   [ D \ (U \P**T * B) ]
        //
        i = n;
        while (i >= 1) {
            if (ipiv[i - 1] > 0) {
                Rscal(nrhs, one / a[(i - 1) + (i - 1) * lda], b[(i - 1)], ldb);
            } else if (i > 1) {
                akm1k = e[i - 1];
                akm1 = a[((i - 1) - 1) + ((i - 1) - 1) * lda] / akm1k;
                ak = a[(i - 1) + (i - 1) * lda] / akm1k;
                denom = akm1 * ak - one;
                for (j = 1; j <= nrhs; j = j + 1) {
                    bkm1 = b[((i - 1) - 1) + (j - 1) * ldb] / akm1k;
                    bk = b[(i - 1) + (j - 1) * ldb] / akm1k;
                    b[((i - 1) - 1) + (j - 1) * ldb] = (ak * bkm1 - bk) / denom;
                    b[(i - 1) + (j - 1) * ldb] = (akm1 * bk - bkm1) / denom;
                }
                i = i - 1;
            }
            i = i - 1;
        }
        //
        //        Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
        //
        Rtrsm("L", "U", "T", "U", n, nrhs, one, a, lda, b, ldb);
        //
        //        P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
        //
        //        Interchange rows K and IPIV(K) of matrix B in reverse order
        //        from the formation order of IPIV(I) vector for Upper case.
        //
        //        (We can do the simple loop over IPIV with increment 1,
        //        since the ABS value of IPIV(I) represents the row index
        //        of the INTEGERerchange with row i in both 1x1 and 2x2 pivot cases)
        //
        for (k = 1; k <= n; k = k + 1) {
            kp = abs(ipiv[k - 1]);
            if (kp != k) {
                Rswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
            }
        }
        //
    } else {
        //
        //        Begin Lower
        //
        //        Solve A*X = B, where A = L*D*L**T.
        //
        //        P**T * B
        //        Interchange rows K and IPIV(K) of matrix B in the same order
        //        that the formation order of IPIV(I) vector for Lower case.
        //
        //        (We can do the simple loop over IPIV with increment 1,
        //        since the ABS value of IPIV(I) represents the row index
        //        of the INTEGERerchange with row i in both 1x1 and 2x2 pivot cases)
        //
        for (k = 1; k <= n; k = k + 1) {
            kp = abs(ipiv[k - 1]);
            if (kp != k) {
                Rswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
            }
        }
        //
        //        Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
        //
        Rtrsm("L", "L", "N", "U", n, nrhs, one, a, lda, b, ldb);
        //
        //        Compute D \ B -> B   [ D \ (L \P**T * B) ]
        //
        i = 1;
        while (i <= n) {
            if (ipiv[i - 1] > 0) {
                Rscal(nrhs, one / a[(i - 1) + (i - 1) * lda], b[(i - 1)], ldb);
            } else if (i < n) {
                akm1k = e[i - 1];
                akm1 = a[(i - 1) + (i - 1) * lda] / akm1k;
                ak = a[((i + 1) - 1) + ((i + 1) - 1) * lda] / akm1k;
                denom = akm1 * ak - one;
                for (j = 1; j <= nrhs; j = j + 1) {
                    bkm1 = b[(i - 1) + (j - 1) * ldb] / akm1k;
                    bk = b[((i + 1) - 1) + (j - 1) * ldb] / akm1k;
                    b[(i - 1) + (j - 1) * ldb] = (ak * bkm1 - bk) / denom;
                    b[((i + 1) - 1) + (j - 1) * ldb] = (akm1 * bk - bkm1) / denom;
                }
                i++;
            }
            i++;
        }
        //
        //        Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
        //
        Rtrsm("L", "L", "T", "U", n, nrhs, one, a, lda, b, ldb);
        //
        //        P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
        //
        //        Interchange rows K and IPIV(K) of matrix B in reverse order
        //        from the formation order of IPIV(I) vector for Lower case.
        //
        //        (We can do the simple loop over IPIV with decrement -1,
        //        since the ABS value of IPIV(I) represents the row index
        //        of the INTEGERerchange with row i in both 1x1 and 2x2 pivot cases)
        //
        for (k = n; k >= 1; k = k - 1) {
            kp = abs(ipiv[k - 1]);
            if (kp != k) {
                Rswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
            }
        }
        //
        //        END Lower
        //
    }
    //
    //     End of Rsytrs_3
    //
}
