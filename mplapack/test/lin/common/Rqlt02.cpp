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

void Rqlt02(common &cmn, INTEGER const m, INTEGER const n, INTEGER const k, REAL *a, REAL *af, REAL *q, REAL *l, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (m == 0 || n == 0 || k == 0) {
        result[1 - 1] = zero;
        result[2 - 1] = zero;
        return;
    }
    //
    REAL eps = Rlamch("Epsilon");
    //
    //     Copy the last k columns of the factorization to the array Q
    //
    const REAL rogue = -1.0e+10;
    dlaset("Full", m, n, rogue, rogue, q, lda);
    if (k < m) {
        dlacpy("Full", m - k, k, af[((n - k + 1) - 1) * ldaf], lda, &q[((n - k + 1) - 1) * ldq], lda);
    }
    if (k > 1) {
        dlacpy("Upper", k - 1, k - 1, af[((m - k + 1) - 1) + ((n - k + 2) - 1) * ldaf], lda, &q[((m - k + 1) - 1) + ((n - k + 2) - 1) * ldq], lda);
    }
    //
    //     Generate the last n columns of the matrix Q
    //
    cmn.srnamt = "DORGQL";
    INTEGER info = 0;
    dorgql(m, n, k, q, lda, &tau[(n - k + 1) - 1], work, lwork, info);
    //
    //     Copy L(m-n+1:m,n-k+1:n)
    //
    dlaset("Full", n, k, zero, zero, l[((m - n + 1) - 1) + ((n - k + 1) - 1) * ldl], lda);
    dlacpy("Lower", k, k, af[((m - k + 1) - 1) + ((n - k + 1) - 1) * ldaf], lda, l[((m - k + 1) - 1) + ((n - k + 1) - 1) * ldl], lda);
    //
    //     Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)
    //
    const REAL one = 1.0;
    Rgemm("Transpose", "No transpose", n, k, m, -one, q, lda, &a[((n - k + 1) - 1) * lda], lda, one, l[((m - n + 1) - 1) + ((n - k + 1) - 1) * ldl], lda);
    //
    //     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
    //
    REAL anorm = dlange("1", m, k, &a[((n - k + 1) - 1) * lda], lda, rwork);
    REAL resid = dlange("1", n, k, l[((m - n + 1) - 1) + ((n - k + 1) - 1) * ldl], lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max((INTEGER)1, m)).real()) / anorm) / eps;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute I - Q'*Q
    //
    dlaset("Full", n, n, zero, one, l, lda);
    Rsyrk("Upper", "Transpose", n, m, -one, q, lda, one, l, lda);
    //
    //     Compute norm( I - Q'*Q ) / ( M * EPS ) .
    //
    resid = dlansy("1", "Upper", n, l, lda, rwork);
    //
    result[2 - 1] = (resid / (max((INTEGER)1, m)).real()) / eps;
    //
    //     End of Rqlt02
    //
}
