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

void Rpot03(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *ainv, INTEGER const ldainv, REAL *work, INTEGER const ldwork, REAL *rwork, REAL &rcond, REAL &resid) {
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
    REAL anorm = Rlansy("1", uplo, n, a, lda, rwork);
    REAL ainvnm = Rlansy("1", uplo, n, ainv, ldainv, rwork);
    if (anorm <= zero || ainvnm <= zero) {
        rcond = zero;
        resid = one / eps;
        return;
    }
    rcond = (one / anorm) / ainvnm;
    //
    //     Expand AINV into a full matrix and call Rsymm to multiply
    //     AINV on the left by A.
    //
    INTEGER j = 0;
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                ainv[(j - 1) + (i - 1) * ldainv] = ainv[(i - 1) + (j - 1) * ldainv];
            }
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            for (i = j + 1; i <= n; i = i + 1) {
                ainv[(j - 1) + (i - 1) * ldainv] = ainv[(i - 1) + (j - 1) * ldainv];
            }
        }
    }
    Rsymm("Left", uplo, n, n, -one, a, lda, ainv, ldainv, zero, work, ldwork);
    //
    //     Add the identity matrix to WORK .
    //
    for (i = 1; i <= n; i = i + 1) {
        work[(i - 1) + (i - 1) * ldwork] += one;
    }
    //
    //     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
    //
    resid = Rlange("1", n, n, work, ldwork, rwork);
    //
    resid = ((resid * rcond) / eps) / n.real();
    //
    //     End of Rpot03
    //
}
