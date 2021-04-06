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

void Claset(const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX const beta, COMPLEX *a, INTEGER const lda) {
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
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER j = 0;
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        //
        //        Set the diagonal to BETA and the strictly upper triangular
        //        part of the array to ALPHA.
        //
        for (j = 2; j <= n; j = j + 1) {
            for (i = 1; i <= min(j - 1, m); i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = alpha;
            }
        }
        for (i = 1; i <= min(n, m); i = i + 1) {
            a[(i - 1) + (i - 1) * lda] = beta;
        }
        //
    } else if (Mlsame(uplo, "L")) {
        //
        //        Set the diagonal to BETA and the strictly lower triangular
        //        part of the array to ALPHA.
        //
        for (j = 1; j <= min(m, n); j = j + 1) {
            for (i = j + 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = alpha;
            }
        }
        for (i = 1; i <= min(n, m); i = i + 1) {
            a[(i - 1) + (i - 1) * lda] = beta;
        }
        //
    } else {
        //
        //        Set the array to BETA on the diagonal and ALPHA on the
        //        offdiagonal.
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = alpha;
            }
        }
        for (i = 1; i <= min(m, n); i = i + 1) {
            a[(i - 1) + (i - 1) * lda] = beta;
        }
    }
    //
    //     End of Claset
    //
}
