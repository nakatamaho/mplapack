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

void Chpt01(const char *uplo, INTEGER const n, COMPLEX *a, COMPLEX *afac, INTEGER *ipiv, COMPLEX *c, INTEGER const ldc, REAL *rwork, REAL &resid) {
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
    REAL anorm = Clanhp("1", uplo, n, a, rwork);
    //
    //     Check the imaginary parts of the diagonal elements and return with
    //     an error code if any are nonzero.
    //
    INTEGER jc = 1;
    INTEGER j = 0;
    const REAL one = 1.0;
    if (Mlsame(uplo, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            if (afac[jc - 1].imag() != zero) {
                resid = one / eps;
                return;
            }
            jc += j + 1;
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            if (afac[jc - 1].imag() != zero) {
                resid = one / eps;
                return;
            }
            jc += n - j + 1;
        }
    }
    //
    //     Initialize C to the identity matrix.
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    Claset("Full", n, n, czero, cone, c, ldc);
    //
    //     Call Clavhp to form the product D * U' (or D * L' ).
    //
    INTEGER info = 0;
    Clavhp(uplo, "Conjugate", "Non-unit", n, n, afac, ipiv, c, ldc, info);
    //
    //     Call Clavhp again to multiply by U ( or L ).
    //
    Clavhp(uplo, "No transpose", "Unit", n, n, afac, ipiv, c, ldc, info);
    //
    //     Compute the difference  C - A .
    //
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        jc = 0;
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(jc + i) - 1];
            }
            c[(j - 1) + (j - 1) * ldc] = c[(j - 1) + (j - 1) * ldc] - a[(jc + j) - 1].real();
            jc += j;
        }
    } else {
        jc = 1;
        for (j = 1; j <= n; j = j + 1) {
            c[(j - 1) + (j - 1) * ldc] = c[(j - 1) + (j - 1) * ldc] - a[jc - 1].real();
            for (i = j + 1; i <= n; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(jc + i - j) - 1];
            }
            jc += n - j + 1;
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
    //     End of Chpt01
    //
}
