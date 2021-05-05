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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Chet01_aa(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *afac, INTEGER const ldafac, INTEGER *ipiv, COMPLEX *c, INTEGER const ldc, REAL *rwork, REAL &resid) {
    a([lda * star]);
    afac([ldafac * star]);
    c([ldc * star]);
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
    //     Initialize C to the tridiagonal matrix T.
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Claset("Full", n, n, czero, czero, c, ldc);
    Clacpy("F", 1, n, afac[(1 - 1)], ldafac + 1, &c[(1 - 1)], ldc + 1);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER j = 0;
    INTEGER i = 0;
    if (n > 1) {
        if (Mlsame(uplo, "U")) {
            Clacpy("F", 1, n - 1, afac[(2 - 1) * ldafac], ldafac + 1, &c[(2 - 1) * ldc], ldc + 1);
            Clacpy("F", 1, n - 1, afac[(2 - 1) * ldafac], ldafac + 1, &c[(2 - 1)], ldc + 1);
            Clacgv(n - 1, &c[(2 - 1)], ldc + 1);
        } else {
            Clacpy("F", 1, n - 1, afac[(2 - 1)], ldafac + 1, &c[(2 - 1) * ldc], ldc + 1);
            Clacpy("F", 1, n - 1, afac[(2 - 1)], ldafac + 1, &c[(2 - 1)], ldc + 1);
            Clacgv(n - 1, &c[(2 - 1) * ldc], ldc + 1);
        }
        //
        //        Call Ctrmm to form the product U' * D (or L * D ).
        //
        if (Mlsame(uplo, "U")) {
            Ctrmm("Left", uplo, "Conjugate transpose", "Unit", n - 1, n, cone, afac[(2 - 1) * ldafac], ldafac, &c[(2 - 1)], ldc);
        } else {
            Ctrmm("Left", uplo, "No transpose", "Unit", n - 1, n, cone, afac[(2 - 1)], ldafac, &c[(2 - 1)], ldc);
        }
        //
        //        Call Ctrmm again to multiply by U (or L ).
        //
        if (Mlsame(uplo, "U")) {
            Ctrmm("Right", uplo, "No transpose", "Unit", n, n - 1, cone, afac[(2 - 1) * ldafac], ldafac, &c[(2 - 1) * ldc], ldc);
        } else {
            Ctrmm("Right", uplo, "Conjugate transpose", "Unit", n, n - 1, cone, afac[(2 - 1)], ldafac, &c[(2 - 1) * ldc], ldc);
        }
        //
        //        Apply hermitian pivots
        //
        for (j = n; j >= 1; j = j - 1) {
            i = ipiv[j - 1];
            if (i != j) {
                Cswap(n, &c[(j - 1)], ldc, &c[(i - 1)], ldc);
            }
        }
        for (j = n; j >= 1; j = j - 1) {
            i = ipiv[j - 1];
            if (i != j) {
                Cswap(n, &c[(j - 1) * ldc], 1, &c[(i - 1) * ldc], 1);
            }
        }
    }
    //
    //     Compute the difference  C - A .
    //
    if (Mlsame(uplo, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(i - 1) + (j - 1) * lda];
            }
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            for (i = j; i <= n; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - a[(i - 1) + (j - 1) * lda];
            }
        }
    }
    //
    //     Compute norm( C - A ) / ( N * norm(A) * EPS )
    //
    resid = Clanhe("1", uplo, n, c, ldc, rwork);
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
    //     End of Chet01_aa
    //
}
