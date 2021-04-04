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

void Rpbstf(const char *uplo, INTEGER const &n, INTEGER const &kd, REAL *ab, INTEGER const &ldab, INTEGER &info) {
    bool upper = false;
    INTEGER kld = 0;
    INTEGER m = 0;
    INTEGER j = 0;
    REAL ajj = 0.0;
    const REAL zero = 0.0;
    INTEGER km = 0;
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
        Mxerbla("Rpbstf", -info);
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
    //     Set the splitting poINTEGER m.
    //
    m = (n + kd) / 2;
    //
    if (upper) {
        //
        //        Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).
        //
        for (j = n; j >= m + 1; j = j - 1) {
            //
            //           Compute s(j,j) and test for non-positive-definiteness.
            //
            ajj = ab[((kd + 1) - 1) + (j - 1) * ldab];
            if (ajj <= zero) {
                goto statement_50;
            }
            ajj = sqrt(ajj);
            ab[((kd + 1) - 1) + (j - 1) * ldab] = ajj;
            km = min(j - 1, kd);
            //
            //           Compute elements j-km:j-1 of the j-th column and update the
            //           the leading submatrix within the band.
            //
            Rscal(km, one / ajj, ab[((kd + 1 - km) - 1) + (j - 1) * ldab], 1);
            Rsyr("Upper", km, -one, ab[((kd + 1 - km) - 1) + (j - 1) * ldab], 1, ab[((kd + 1) - 1) + ((j - km) - 1) * ldab], kld);
        }
        //
        //        Factorize the updated submatrix A(1:m,1:m) as U**T*U.
        //
        for (j = 1; j <= m; j = j + 1) {
            //
            //           Compute s(j,j) and test for non-positive-definiteness.
            //
            ajj = ab[((kd + 1) - 1) + (j - 1) * ldab];
            if (ajj <= zero) {
                goto statement_50;
            }
            ajj = sqrt(ajj);
            ab[((kd + 1) - 1) + (j - 1) * ldab] = ajj;
            km = min(kd, m - j);
            //
            //           Compute elements j+1:j+km of the j-th row and update the
            //           trailing submatrix within the band.
            //
            if (km > 0) {
                Rscal(km, one / ajj, ab[(kd - 1) + ((j + 1) - 1) * ldab], kld);
                Rsyr("Upper", km, -one, ab[(kd - 1) + ((j + 1) - 1) * ldab], kld, ab[((kd + 1) - 1) + ((j + 1) - 1) * ldab], kld);
            }
        }
    } else {
        //
        //        Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).
        //
        for (j = n; j >= m + 1; j = j - 1) {
            //
            //           Compute s(j,j) and test for non-positive-definiteness.
            //
            ajj = ab[(j - 1) * ldab];
            if (ajj <= zero) {
                goto statement_50;
            }
            ajj = sqrt(ajj);
            ab[(j - 1) * ldab] = ajj;
            km = min(j - 1, kd);
            //
            //           Compute elements j-km:j-1 of the j-th row and update the
            //           trailing submatrix within the band.
            //
            Rscal(km, one / ajj, ab[((km + 1) - 1) + ((j - km) - 1) * ldab], kld);
            Rsyr("Lower", km, -one, ab[((km + 1) - 1) + ((j - km) - 1) * ldab], kld, ab[((j - km) - 1) * ldab], kld);
        }
        //
        //        Factorize the updated submatrix A(1:m,1:m) as U**T*U.
        //
        for (j = 1; j <= m; j = j + 1) {
            //
            //           Compute s(j,j) and test for non-positive-definiteness.
            //
            ajj = ab[(j - 1) * ldab];
            if (ajj <= zero) {
                goto statement_50;
            }
            ajj = sqrt(ajj);
            ab[(j - 1) * ldab] = ajj;
            km = min(kd, m - j);
            //
            //           Compute elements j+1:j+km of the j-th column and update the
            //           trailing submatrix within the band.
            //
            if (km > 0) {
                Rscal(km, one / ajj, ab[(2 - 1) + (j - 1) * ldab], 1);
                Rsyr("Lower", km, -one, ab[(2 - 1) + (j - 1) * ldab], 1, ab[((j + 1) - 1) * ldab], kld);
            }
        }
    }
    return;
//
statement_50:
    info = j;
    //
    //     End of Rpbstf
    //
}
