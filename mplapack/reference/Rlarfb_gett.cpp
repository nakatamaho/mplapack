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

void Rlarfb_gett(const char *ident, INTEGER const m, INTEGER const n, INTEGER const k, REAL *t, INTEGER const ldt, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *work, INTEGER const ldwork) {
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. EXTERNAL FUNCTIONS ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (m < 0 || n <= 0 || k == 0 || k > n) {
        return;
    }
    //
    bool lnotident = !Mlsame(ident, "I");
    //
    //     ------------------------------------------------------------------
    //
    //     First Step. Computation of the Column Block 2:
    //
    //        ( A2 ) := H * ( A2 )
    //        ( B2 )        ( B2 )
    //
    //     ------------------------------------------------------------------
    //
    INTEGER j = 0;
    const REAL one = 1.0;
    INTEGER i = 0;
    if (n > k) {
        //
        //        col2_(1) Compute W2: = A2. Therefore, copy A2 = A(1:K, K+1:N)
        //        INTEGERo W2=WORK(1:K, 1:N-K) column-by-column.
        //
        for (j = 1; j <= n - k; j = j + 1) {
            Rcopy(k, &a[((k + j) - 1) * lda], 1, &work[(j - 1) * ldwork], 1);
        }
        //
        if (lnotident) {
            //
            //           col2_(2) Compute W2: = (V1**T) * W2 = (A1**T) * W2,
            //           V1 is not an identy matrix, but unit lower-triangular
            //           V1 stored in A1 (diagonal ones are not stored).
            //
            Rtrmm("L", "L", "T", "U", k, n - k, one, a, lda, work, ldwork);
        }
        //
        //        col2_(3) Compute W2: = W2 + (V2**T) * B2 = W2 + (B1**T) * B2
        //        V2 stored in B1.
        //
        if (m > 0) {
            Rgemm("T", "N", k, n - k, m, one, b, ldb, &b[((k + 1) - 1) * ldb], ldb, one, work, ldwork);
        }
        //
        //        col2_(4) Compute W2: = T * W2,
        //        T is upper-triangular.
        //
        Rtrmm("L", "U", "N", "N", k, n - k, one, t, ldt, work, ldwork);
        //
        //        col2_(5) Compute B2: = B2 - V2 * W2 = B2 - B1 * W2,
        //        V2 stored in B1.
        //
        if (m > 0) {
            Rgemm("N", "N", m, n - k, k, -one, b, ldb, work, ldwork, one, &b[((k + 1) - 1) * ldb], ldb);
        }
        //
        if (lnotident) {
            //
            //           col2_(6) Compute W2: = V1 * W2 = A1 * W2,
            //           V1 is not an identity matrix, but unit lower-triangular,
            //           V1 stored in A1 (diagonal ones are not stored).
            //
            Rtrmm("L", "L", "N", "U", k, n - k, one, a, lda, work, ldwork);
        }
        //
        //        col2_(7) Compute A2: = A2 - W2 =
        //                             = A(1:K, K+1:N-K) - WORK(1:K, 1:N-K),
        //        column-by-column.
        //
        for (j = 1; j <= n - k; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                a[(i - 1) + ((k + j) - 1) * lda] = a[(i - 1) + ((k + j) - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
    }
    //
    //     ------------------------------------------------------------------
    //
    //     Second Step. Computation of the Column Block 1:
    //
    //        ( A1 ) := H * ( A1 )
    //        ( B1 )        (  0 )
    //
    //     ------------------------------------------------------------------
    //
    //     col1_(1) Compute W1: = A1. Copy the upper-triangular
    //     A1 = A(1:K, 1:K) INTEGERo the upper-triangular
    //     W1 = WORK(1:K, 1:K) column-by-column.
    //
    for (j = 1; j <= k; j = j + 1) {
        Rcopy(j, &a[(j - 1) * lda], 1, &work[(j - 1) * ldwork], 1);
    }
    //
    //     Set the subdiagonal elements of W1 to zero column-by-column.
    //
    const REAL zero = 0.0;
    for (j = 1; j <= k - 1; j = j + 1) {
        for (i = j + 1; i <= k; i = i + 1) {
            work[(i - 1) + (j - 1) * ldwork] = zero;
        }
    }
    //
    if (lnotident) {
        //
        //        col1_(2) Compute W1: = (V1**T) * W1 = (A1**T) * W1,
        //        V1 is not an identity matrix, but unit lower-triangular
        //        V1 stored in A1 (diagonal ones are not stored),
        //        W1 is upper-triangular with zeroes below the diagonal.
        //
        Rtrmm("L", "L", "T", "U", k, k, one, a, lda, work, ldwork);
    }
    //
    //     col1_(3) Compute W1: = T * W1,
    //     T is upper-triangular,
    //     W1 is upper-triangular with zeroes below the diagonal.
    //
    Rtrmm("L", "U", "N", "N", k, k, one, t, ldt, work, ldwork);
    //
    //     col1_(4) Compute B1: = - V2 * W1 = - B1 * W1,
    //     V2 = B1, W1 is upper-triangular with zeroes below the diagonal.
    //
    if (m > 0) {
        Rtrmm("R", "U", "N", "N", m, k, -one, work, ldwork, b, ldb);
    }
    //
    if (lnotident) {
        //
        //        col1_(5) Compute W1: = V1 * W1 = A1 * W1,
        //        V1 is not an identity matrix, but unit lower-triangular
        //        V1 stored in A1 (diagonal ones are not stored),
        //        W1 is upper-triangular on input with zeroes below the diagonal,
        //        and square on output.
        //
        Rtrmm("L", "L", "N", "U", k, k, one, a, lda, work, ldwork);
        //
        //        col1_(6) Compute A1: = A1 - W1 = A(1:K, 1:K) - WORK(1:K, 1:K)
        //        column-by-column. A1 is upper-triangular on input.
        //        If IDENT, A1 is square on output, and W1 is square,
        //        if NOT IDENT, A1 is upper-triangular on output,
        //        W1 is upper-triangular.
        //
        //        col1_(6)_a Compute elements of A1 below the diagonal.
        //
        for (j = 1; j <= k - 1; j = j + 1) {
            for (i = j + 1; i <= k; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = -work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
    }
    //
    //     col1_(6)_b Compute elements of A1 on and above the diagonal.
    //
    for (j = 1; j <= k; j = j + 1) {
        for (i = 1; i <= j; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - work[(i - 1) + (j - 1) * ldwork];
        }
    }
    //
    //     End of Rlarfb_gett
    //
}
