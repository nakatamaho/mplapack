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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

void Chet01_rook(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *afac, INTEGER const ldafac, INTEGER *ipiv, COMPLEX *c, INTEGER const ldc, REAL *rwork, REAL &resid) {
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
    REAL anorm = Clanhe("1", uplo, n, a, lda, rwork);
    //
    //     Check the imaginary parts of the diagonal elements and return with
    //     an error code if any are nonzero.
    //
    INTEGER j = 0;
    const REAL one = 1.0;
    for (j = 1; j <= n; j = j + 1) {
        if (afac[(j - 1) + (j - 1) * ldafac].imag() != zero) {
            resid = one / eps;
            return;
        }
    }
    //
    //     Initialize C to the identity matrix.
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    Claset("Full", n, n, czero, cone, c, ldc);
    //
    //     Call Clavhe_rook to form the product D * U' (or D * L' ).
    //
    INTEGER info = 0;
    Clavhe_rook(uplo, "Conjugate", "Non-unit", n, n, afac, ldafac, ipiv, c, ldc, info);
    //
    //     Call Clavhe_rook again to multiply by U (or L ).
    //
    Clavhe_rook(uplo, "No transpose", "Unit", n, n, afac, ldafac, ipiv, c, ldc, info);
    //
    //     Compute the difference  C - A .
    //
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(i - 1) + (j - 1) * lda];
            }
            c[(j - 1) + (j - 1) * ldc] = c[(j - 1) + (j - 1) * ldc] - a[(j - 1) + (j - 1) * lda].real();
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            c[(j - 1) + (j - 1) * ldc] = c[(j - 1) + (j - 1) * ldc] - a[(j - 1) + (j - 1) * lda].real();
            for (i = j + 1; i <= n; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(i - 1) + (j - 1) * lda];
            }
        }
    }
    //
    //     Compute norm( C - A ) / ( N * norm(A) * EPS )
    //
    resid = Clanhe("1", uplo, n, c, ldc, rwork);
    //
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        resid = ((resid / n.real()) / anorm) / eps;
    }
    //
    //     End of Chet01_rook
    //
}
