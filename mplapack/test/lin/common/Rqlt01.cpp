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

#include <mplapack_debug.h>

void Rqlt01(INTEGER const m, INTEGER const n, REAL *a, REAL *af, REAL *q, REAL *l, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    INTEGER ldaf = lda;
    INTEGER ldq = lda;
    INTEGER ldl = lda;
    INTEGER minmn = min(m, n);
    REAL eps = Rlamch("Epsilon");
    //
    //     Copy the matrix A to the array AF.
    //
    Rlacpy("Full", m, n, a, lda, af, lda);
    //
    //     Factorize the matrix A in the array AF.
    //
    INTEGER info = 0;
    Rgeqlf(m, n, af, lda, tau, work, lwork, info);
    //
    //     Copy details of Q
    //
    const REAL rogue = -1.0e+10;
    Rlaset("Full", m, m, rogue, rogue, q, lda);
    if (m >= n) {
        if (n < m && n > 0) {
            Rlacpy("Full", m - n, n, af, lda, &q[((m - n + 1) - 1) * ldq], lda);
        }
        if (n > 1) {
            Rlacpy("Upper", n - 1, n - 1, &af[((m - n + 1) - 1) + (2 - 1) * ldaf], lda, &q[((m - n + 1) - 1) + ((m - n + 2) - 1) * ldq], lda);
        }
    } else {
        if (m > 1) {
            Rlacpy("Upper", m - 1, m - 1, &af[((n - m + 2) - 1) * ldaf], lda, &q[(2 - 1) * ldq], lda);
        }
    }
    //
    //     Generate the m-by-m matrix Q
    //
    strncpy(srnamt, "Rorgql", srnamt_len);
    Rorgql(m, m, minmn, q, lda, tau, work, lwork, info);
    //
    //     Copy L
    //
    const REAL zero = 0.0;
    Rlaset("Full", m, n, zero, zero, l, lda);
    if (m >= n) {
        if (n > 0) {
            Rlacpy("Lower", n, n, &af[((m - n + 1) - 1)], lda, &l[((m - n + 1) - 1)], lda);
        }
    } else {
        if (n > m && m > 0) {
            Rlacpy("Full", m, n - m, af, lda, l, lda);
        }
        if (m > 0) {
            Rlacpy("Lower", m, m, &af[((n - m + 1) - 1) * ldaf], lda, &l[((n - m + 1) - 1) * ldl], lda);
        }
    }
    //
    //     Compute L - Q'*A
    //
    const REAL one = 1.0;
    Rgemm("Transpose", "No transpose", m, n, m, -one, q, lda, a, lda, one, l, lda);
    //
    //     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
    //
    REAL anorm = Rlange("1", m, n, a, lda, rwork);
    REAL resid = Rlange("1", m, n, l, lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / castREAL(max((INTEGER)1, m))) / anorm) / eps;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute I - Q'*Q
    //
    Rlaset("Full", m, m, zero, one, l, lda);
    Rsyrk("Upper", "Transpose", m, m, -one, q, lda, one, l, lda);
    //
    //     Compute norm( I - Q'*Q ) / ( M * EPS ) .
    //
    resid = Rlansy("1", "Upper", m, l, lda, rwork);
    //
    result[2 - 1] = (resid / castREAL(max((INTEGER)1, m))) / eps;
    //
    //     End of Rqlt01
    //
}
