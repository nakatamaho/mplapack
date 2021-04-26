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

void Cspt01(const char *uplo, INTEGER const n, COMPLEX *a, COMPLEX *afac, INTEGER *ipiv, COMPLEX *c, INTEGER const ldc, REAL *rwork, REAL &resid) {
    //
    //  -- LAPACK test routine --
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
    //     Quick exit if N = 0.
    //
    const REAL zero = 0.0;
    if (n <= 0) {
        resid = zero;
        return;
    }
    //
    //     Determine EPS and the norm of A.
    //
    REAL eps = Rlamch("Epsilon");
    REAL anorm = zlansp("1", uplo, n, a, rwork);
    //
    //     Initialize C to the identity matrix.
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    zlaset("Full", n, n, czero, cone, c, ldc);
    //
    //     Call Clavsp to form the product D * U' (or D * L' ).
    //
    INTEGER info = 0;
    Clavsp(uplo, "Transpose", "Non-unit", n, n, afac, ipiv, c, ldc, info);
    //
    //     Call Clavsp again to multiply by U ( or L ).
    //
    Clavsp(uplo, "No transpose", "Unit", n, n, afac, ipiv, c, ldc, info);
    //
    //     Compute the difference  C - A .
    //
    INTEGER jc = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        jc = 0;
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(jc + i) - 1];
            }
            jc += j;
        }
    } else {
        jc = 1;
        for (j = 1; j <= n; j = j + 1) {
            for (i = j; i <= n; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(jc + i - j) - 1];
            }
            jc += n - j + 1;
        }
    }
    //
    //     Compute norm( C - A ) / ( N * norm(A) * EPS )
    //
    resid = zlansy("1", uplo, n, c, ldc, rwork);
    //
    const REAL one = 1.0;
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        resid = ((resid / n.real()) / anorm) / eps;
    }
    //
    //     End of Cspt01
    //
}
