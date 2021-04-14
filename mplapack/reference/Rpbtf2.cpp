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

void Rpbtf2(const char *uplo, INTEGER const n, INTEGER const kd, REAL *ab, INTEGER const ldab, INTEGER &info) {
    bool upper = false;
    INTEGER kld = 0;
    INTEGER j = 0;
    REAL ajj = 0.0;
    const REAL zero = 0.0;
    INTEGER kn = 0;
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
    } else if (kd < 0) {
        info = -3;
    } else if (ldab < kd + 1) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Rpbtf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    kld = max((INTEGER)1, ldab - 1);
    //
    if (upper) {
        //
        //        Compute the Cholesky factorization A = U**T*U.
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //           Compute U(J,J) and test for non-positive-definiteness.
            //
            ajj = ab[((kd + 1) - 1) + (j - 1) * ldab];
            if (ajj <= zero) {
                goto statement_30;
            }
            ajj = sqrt(ajj);
            ab[((kd + 1) - 1) + (j - 1) * ldab] = ajj;
            //
            //           Compute elements J+1:J+KN of row J and update the
            //           trailing submatrix within the band.
            //
            kn = min(kd, n - j);
            if (kn > 0) {
                Rscal(kn, one / ajj, ab[(kd - 1) + ((j + 1) - 1) * ldab], kld);
                Rsyr("Upper", kn, -one, ab[(kd - 1) + ((j + 1) - 1) * ldab], kld, ab[((kd + 1) - 1) + ((j + 1) - 1) * ldab], kld);
            }
        }
    } else {
        //
        //        Compute the Cholesky factorization A = L*L**T.
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //           Compute L(J,J) and test for non-positive-definiteness.
            //
            ajj = ab[(j - 1) * ldab];
            if (ajj <= zero) {
                goto statement_30;
            }
            ajj = sqrt(ajj);
            ab[(j - 1) * ldab] = ajj;
            //
            //           Compute elements J+1:J+KN of column J and update the
            //           trailing submatrix within the band.
            //
            kn = min(kd, n - j);
            if (kn > 0) {
                Rscal(kn, one / ajj, ab[(2 - 1) + (j - 1) * ldab], 1);
                Rsyr("Lower", kn, -one, ab[(2 - 1) + (j - 1) * ldab], 1, ab[((j + 1) - 1) * ldab], kld);
            }
        }
    }
    return;
//
statement_30:
    info = j;
    //
    //     End of Rpbtf2
    //
}
