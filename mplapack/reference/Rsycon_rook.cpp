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

void Rsycon_rook(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, INTEGER *ipiv, REAL const &anorm, REAL &rcond, REAL *work, INTEGER *iwork, INTEGER &info) {
    bool upper = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER i = 0;
    INTEGER kase = 0;
    REAL ainvnm = 0.0;
    arr_1d<3, INTEGER> isave(fill0);
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
    //     Test the input parameters.
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (anorm < zero) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rsycon_rook", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    rcond = zero;
    if (n == 0) {
        rcond = one;
        return;
    } else if (anorm <= zero) {
        return;
    }
    //
    //     Check that the diagonal matrix D is nonsingular.
    //
    if (upper) {
        //
        //        Upper triangular storage: examine D from bottom to top
        //
        for (i = n; i >= 1; i = i - 1) {
            if (ipiv[i - 1] > 0 && a[(i - 1) + (i - 1) * lda] == zero) {
                return;
            }
        }
    } else {
        //
        //        Lower triangular storage: examine D from top to bottom.
        //
        for (i = 1; i <= n; i = i + 1) {
            if (ipiv[i - 1] > 0 && a[(i - 1) + (i - 1) * lda] == zero) {
                return;
            }
        }
    }
    //
    //     Estimate the 1-norm of the inverse.
    //
    kase = 0;
statement_30:
    Rlacn2(n, work[(n + 1) - 1], work, iwork, ainvnm, kase, isave);
    if (kase != 0) {
        //
        //        Multiply by inv(L*D*L**T) or inv(U*D*U**T).
        //
        Rsytrs_rook(uplo, n, 1, a, lda, ipiv, work, n, info);
        goto statement_30;
    }
    //
    //     Compute the estimate of the reciprocal condition number.
    //
    if (ainvnm != zero) {
        rcond = (one / ainvnm) / anorm;
    }
    //
    //     End of Rsycon_rook
    //
}
