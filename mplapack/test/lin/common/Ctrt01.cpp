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

void Ctrt01(const char *uplo, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *ainv, INTEGER const ldainv, REAL &rcond, REAL *rwork, REAL &resid) {
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
    //     Quick exit if N = 0
    //
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (n <= 0) {
        rcond = one;
        resid = zero;
        return;
    }
    //
    //     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    REAL anorm = Clantr("1", uplo, diag, n, n, a, lda, rwork);
    REAL ainvnm = Clantr("1", uplo, diag, n, n, ainv, ldainv, rwork);
    if (anorm <= zero || ainvnm <= zero) {
        rcond = zero;
        resid = one / eps;
        return;
    }
    rcond = (one / anorm) / ainvnm;
    //
    //     Set the diagonal of AINV to 1 if AINV has unit diagonal.
    //
    INTEGER j = 0;
    if (Mlsame(diag, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            ainv[(j - 1) + (j - 1) * ldainv] = one;
        }
    }
    //
    //     Compute A * AINV, overwriting AINV.
    //
    if (Mlsame(uplo, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            Ctrmv("Upper", "No transpose", diag, j, a, lda, ainv[(j - 1) * ldainv], 1);
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            Ctrmv("Lower", "No transpose", diag, n - j + 1, &a[(j - 1) + (j - 1) * lda], lda, ainv[(j - 1) + (j - 1) * ldainv], 1);
        }
    }
    //
    //     Subtract 1 from each diagonal element to form A*AINV - I.
    //
    for (j = 1; j <= n; j = j + 1) {
        ainv[(j - 1) + (j - 1) * ldainv] = ainv[(j - 1) + (j - 1) * ldainv] - one;
    }
    //
    //     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
    //
    resid = Clantr("1", uplo, "Non-unit", n, n, ainv, ldainv, rwork);
    //
    resid = ((resid * rcond) / n.real()) / eps;
    //
    //     End of Ctrt01
    //
}
