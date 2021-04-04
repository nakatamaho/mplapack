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

REAL Clanhb(const char *norm, const char *uplo, INTEGER const &n, INTEGER const &k, COMPLEX *ab, INTEGER const &ldab, REAL *work) {
    REAL return_value = 0.0;
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
    // =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    REAL value = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    REAL sum = 0.0;
    INTEGER l = 0;
    REAL absa = 0.0;
    arr_1d<2, REAL> ssq(fill0);
    const REAL one = 1.0;
    arr_1d<2, REAL> colssq(fill0);
    if (n == 0) {
        value = zero;
    } else if (Mlsame(norm, "M")) {
        //
        //        Find max(abs(A(i,j))).
        //
        value = zero;
        if (Mlsame(uplo, "U")) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = max(k + 2 - j, 1); i <= k; i = i + 1) {
                    sum = abs(ab[(i - 1) + (j - 1) * ldab]);
                    if (value < sum || Risnan(sum)) {
                        value = sum;
                    }
                }
                sum = abs(ab[((k + 1) - 1) + (j - 1) * ldab].real());
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                sum = abs(ab[(j - 1) * ldab].real());
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
                for (i = 2; i <= min(n + 1 - j, k + 1); i = i + 1) {
                    sum = abs(ab[(i - 1) + (j - 1) * ldab]);
                    if (value < sum || Risnan(sum)) {
                        value = sum;
                    }
                }
            }
        }
    } else if ((Mlsame(norm, "I")) || (Mlsame(norm, "O")) || (norm == "1")) {
        //
        //        Find normI(A) ( = norm1(A), since A is hermitian).
        //
        value = zero;
        if (Mlsame(uplo, "U")) {
            for (j = 1; j <= n; j = j + 1) {
                sum = zero;
                l = k + 1 - j;
                for (i = max((INTEGER)1, j - k); i <= j - 1; i = i + 1) {
                    absa = abs(ab[((l + i) - 1) + (j - 1) * ldab]);
                    sum += absa;
                    work[i - 1] += absa;
                }
                work[j - 1] = sum + abs(ab[((k + 1) - 1) + (j - 1) * ldab].real());
            }
            for (i = 1; i <= n; i = i + 1) {
                sum = work[i - 1];
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        } else {
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = zero;
            }
            for (j = 1; j <= n; j = j + 1) {
                sum = work[j - 1] + abs(ab[(j - 1) * ldab].real());
                l = 1 - j;
                for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                    absa = abs(ab[((l + i) - 1) + (j - 1) * ldab]);
                    sum += absa;
                    work[i - 1] += absa;
                }
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        }
    } else if ((Mlsame(norm, "F")) || (Mlsame(norm, "E"))) {
        //
        //        Find normF(A).
        //        SSQ(1) is scale
        //        SSQ(2) is sum-of-squares
        //        For better accuracy, sum each column separately.
        //
        ssq[1 - 1] = zero;
        ssq[2 - 1] = one;
        //
        //        Sum off-diagonals
        //
        if (k > 0) {
            if (Mlsame(uplo, "U")) {
                for (j = 2; j <= n; j = j + 1) {
                    colssq[1 - 1] = zero;
                    colssq[2 - 1] = one;
                    Classq(min(j - 1, k), ab[((max(k + 2 - j) - 1) + (1) - 1) * ldab], 1, colssq[1 - 1], colssq[2 - 1]);
                    Rcombssq(ssq, colssq);
                }
                l = k + 1;
            } else {
                for (j = 1; j <= n - 1; j = j + 1) {
                    colssq[1 - 1] = zero;
                    colssq[2 - 1] = one;
                    Classq(min(n - j, k), ab[(2 - 1) + (j - 1) * ldab], 1, colssq[1 - 1], colssq[2 - 1]);
                    Rcombssq(ssq, colssq);
                }
                l = 1;
            }
            ssq[2 - 1] = 2 * ssq[2 - 1];
        } else {
            l = 1;
        }
        //
        //        Sum diagonal
        //
        colssq[1 - 1] = zero;
        colssq[2 - 1] = one;
        for (j = 1; j <= n; j = j + 1) {
            if (ab[(l - 1) + (j - 1) * ldab].real() != zero) {
                absa = abs(ab[(l - 1) + (j - 1) * ldab].real());
                if (colssq[1 - 1] < absa) {
                    colssq[2 - 1] = one + colssq[2 - 1] * pow2((colssq[1 - 1] / absa));
                    colssq[1 - 1] = absa;
                } else {
                    colssq[2 - 1] += pow2((absa / colssq[1 - 1]));
                }
            }
        }
        Rcombssq(ssq, colssq);
        value = ssq[1 - 1] * sqrt(ssq[2 - 1]);
    }
    //
    return_value = value;
    return return_value;
    //
    //     End of Clanhb
    //
}
