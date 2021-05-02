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

void Rqrt02(INTEGER const m, INTEGER const n, INTEGER const k, REAL *a, REAL *af, REAL *q, REAL *r, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL eps = Rlamch("Epsilon");
    //
    //     Copy the first k columns of the factorization to the array Q
    //
    const REAL rogue = -1.0e+10;
    Rlaset("Full", m, n, rogue, rogue, q, lda);
    Rlacpy("Lower", m - 1, k, af[(2 - 1)], lda, &q[(2 - 1)], lda);
    //
    //     Generate the first n columns of the matrix Q
    //
    cmn.srnamt = "Rorgqr";
    INTEGER info = 0;
    Rorgqr(m, n, k, q, lda, tau, work, lwork, info);
    //
    //     Copy R(1:n,1:k)
    //
    const REAL zero = 0.0;
    Rlaset("Full", n, k, zero, zero, r, lda);
    Rlacpy("Upper", n, k, af, lda, r, lda);
    //
    //     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
    //
    const REAL one = 1.0;
    Rgemm("Transpose", "No transpose", n, k, m, -one, q, lda, a, lda, one, r, lda);
    //
    //     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
    //
    REAL anorm = Rlange("1", m, k, a, lda, rwork);
    REAL resid = Rlange("1", n, k, r, lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max((INTEGER)1, m)).real()) / anorm) / eps;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute I - Q'*Q
    //
    Rlaset("Full", n, n, zero, one, r, lda);
    Rsyrk("Upper", "Transpose", n, m, -one, q, lda, one, r, lda);
    //
    //     Compute norm( I - Q'*Q ) / ( M * EPS ) .
    //
    resid = Rlansy("1", "Upper", n, r, lda, rwork);
    //
    result[2 - 1] = (resid / (max((INTEGER)1, m)).real()) / eps;
    //
    //     End of Rqrt02
    //
}
