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

void Rrqt01(INTEGER const m, INTEGER const n, REAL *a, REAL *af, REAL *q, REAL *r, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
    // COMMON srnamc
    char &srnamt = cmn.srnamt;
    //
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
    INTEGER minmn = min(m, n);
    REAL eps = Rlamch("Epsilon");
    //
    //     Copy the matrix A to the array AF.
    //
    Rlacpy("Full", m, n, a, lda, af, lda);
    //
    //     Factorize the matrix A in the array AF.
    //
    srnamt = "RgerQF";
    INTEGER info = 0;
    Rgerqf(m, n, af, lda, tau, work, lwork, info);
    //
    //     Copy details of Q
    //
    const REAL rogue = -1.0e+10;
    Rlaset("Full", n, n, rogue, rogue, q, lda);
    if (m <= n) {
        if (m > 0 && m < n) {
            Rlacpy("Full", m, n - m, af, lda, &q[((n - m + 1) - 1)], lda);
        }
        if (m > 1) {
            Rlacpy("Lower", m - 1, m - 1, af[(2 - 1) + ((n - m + 1) - 1) * ldaf], lda, &q[((n - m + 2) - 1) + ((n - m + 1) - 1) * ldq], lda);
        }
    } else {
        if (n > 1) {
            Rlacpy("Lower", n - 1, n - 1, af[((m - n + 2) - 1)], lda, &q[(2 - 1)], lda);
        }
    }
    //
    //     Generate the n-by-n matrix Q
    //
    srnamt = "Rorgrq";
    Rorgrq(n, n, minmn, q, lda, tau, work, lwork, info);
    //
    //     Copy R
    //
    const REAL zero = 0.0;
    Rlaset("Full", m, n, zero, zero, r, lda);
    if (m <= n) {
        if (m > 0) {
            Rlacpy("Upper", m, m, af[((n - m + 1) - 1) * ldaf], lda, r[((n - m + 1) - 1) * ldr], lda);
        }
    } else {
        if (m > n && n > 0) {
            Rlacpy("Full", m - n, n, af, lda, r, lda);
        }
        if (n > 0) {
            Rlacpy("Upper", n, n, af[((m - n + 1) - 1)], lda, r[((m - n + 1) - 1)], lda);
        }
    }
    //
    //     Compute R - A*Q'
    //
    const REAL one = 1.0;
    Rgemm("No transpose", "Transpose", m, n, n, -one, a, lda, q, lda, one, r, lda);
    //
    //     Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .
    //
    REAL anorm = Rlange("1", m, n, a, lda, rwork);
    REAL resid = Rlange("1", m, n, r, lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max((INTEGER)1, n)).real()) / anorm) / eps;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute I - Q*Q'
    //
    Rlaset("Full", n, n, zero, one, r, lda);
    Rsyrk("Upper", "No transpose", n, n, -one, q, lda, one, r, lda);
    //
    //     Compute norm( I - Q*Q' ) / ( N * EPS ) .
    //
    resid = Rlansy("1", "Upper", n, r, lda, rwork);
    //
    result[2 - 1] = (resid / (max((INTEGER)1, n)).real()) / eps;
    //
    //     End of Rrqt01
    //
}
