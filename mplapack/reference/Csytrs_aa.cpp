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

void Csytrs_aa(const char *uplo, INTEGER const &n, INTEGER const &nrhs, COMPLEX *a, INTEGER const &lda, INTEGER *ipiv, COMPLEX *b, INTEGER const &ldb, COMPLEX *work, INTEGER const &lwork, INTEGER &info) {
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
    bool lquery = (lwork == -1);
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
    } else if (lwork < max((INTEGER)1, 3 * n - 2) && !lquery) {
        info = -10;
    }
    INTEGER lwkopt = 0;
    if (info != 0) {
        Mxerbla("Csytrs_aa", -info);
        return;
    } else if (lquery) {
        lwkopt = (3 * n - 2);
        work[1 - 1] = lwkopt;
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
    const COMPLEX one = 1.0;
    if (upper) {
        //
        //        Solve A*X = B, where A = U**T*T*U.
        //
        //        1) Forward substitution with U**T
        //
        if (n > 1) {
            //
            //           Pivot, P**T * B -> B
            //
            for (k = 1; k <= n; k = k + 1) {
                kp = ipiv[k - 1];
                if (kp != k) {
                    Cswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
                }
            }
            //
            //           Compute U**T \ B -> B    [ (U**T \P**T * B) ]
            //
            Ctrsm("L", "U", "T", "U", n - 1, nrhs, one, a[(2 - 1) * lda], lda, b[(2 - 1)], ldb);
        }
        //
        //        2) Solve with triangular matrix T
        //
        //        Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
        //
        Clacpy("F", 1, n, a[(1 - 1)], lda + 1, work[n - 1], 1);
        if (n > 1) {
            Clacpy("F", 1, n - 1, a[(2 - 1) * lda], lda + 1, work[1 - 1], 1);
            Clacpy("F", 1, n - 1, a[(2 - 1) * lda], lda + 1, work[(2 * n) - 1], 1);
        }
        Cgtsv(n, nrhs, work[1 - 1], work[n - 1], work[(2 * n) - 1], b, ldb, info);
        //
        //        3) Backward substitution with U
        //
        if (n > 1) {
            //
            //           Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]
            //
            Ctrsm("L", "U", "N", "U", n - 1, nrhs, one, a[(2 - 1) * lda], lda, b[(2 - 1)], ldb);
            //
            //           Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]
            //
            for (k = n; k >= 1; k = k - 1) {
                kp = ipiv[k - 1];
                if (kp != k) {
                    Cswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
                }
            }
        }
        //
    } else {
        //
        //        Solve A*X = B, where A = L*T*L**T.
        //
        //        1) Forward substitution with L
        //
        if (n > 1) {
            //
            //           Pivot, P**T * B -> B
            //
            for (k = 1; k <= n; k = k + 1) {
                kp = ipiv[k - 1];
                if (kp != k) {
                    Cswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
                }
            }
            //
            //           Compute L \ B -> B    [ (L \P**T * B) ]
            //
            Ctrsm("L", "L", "N", "U", n - 1, nrhs, one, a[(2 - 1)], lda, b[(2 - 1)], ldb);
        }
        //
        //        2) Solve with triangular matrix T
        //
        //        Compute T \ B -> B   [ T \ (L \P**T * B) ]
        //
        Clacpy("F", 1, n, a[(1 - 1)], lda + 1, work[n - 1], 1);
        if (n > 1) {
            Clacpy("F", 1, n - 1, a[(2 - 1)], lda + 1, work[1 - 1], 1);
            Clacpy("F", 1, n - 1, a[(2 - 1)], lda + 1, work[(2 * n) - 1], 1);
        }
        Cgtsv(n, nrhs, work[1 - 1], work[n - 1], work[(2 * n) - 1], b, ldb, info);
        //
        //        3) Backward substitution with L**T
        //
        if (n > 1) {
            //
            //           Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]
            //
            Ctrsm("L", "L", "T", "U", n - 1, nrhs, one, a[(2 - 1)], lda, b[(2 - 1)], ldb);
            //
            //           Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
            //
            for (k = n; k >= 1; k = k - 1) {
                kp = ipiv[k - 1];
                if (kp != k) {
                    Cswap(nrhs, b[(k - 1)], ldb, b[(kp - 1)], ldb);
                }
            }
        }
        //
    }
    //
    //     End of Csytrs_aa
    //
}
