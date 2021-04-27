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

void Clqt02(common &cmn, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, COMPLEX *af, COMPLEX *q, COMPLEX *l, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    //     Copy the first k rows of the factorization to the array Q
    //
    const COMPLEX rogue = COMPLEX(-1.0e+10, -1.0e+10);
    Claset("Full", m, n, rogue, rogue, q, lda);
    Clacpy("Upper", k, n - 1, af[(2 - 1) * ldaf], lda, &q[(2 - 1) * ldq], lda);
    //
    //     Generate the first n columns of the matrix Q
    //
    cmn.srnamt = "Cunglq";
    INTEGER info = 0;
    Cunglq(m, n, k, q, lda, tau, work, lwork, info);
    //
    //     Copy L(1:k,1:m)
    //
    const REAL zero = 0.0;
    Claset("Full", k, m, COMPLEX(zero), COMPLEX(zero), l, lda);
    Clacpy("Lower", k, m, af, lda, l, lda);
    //
    //     Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
    //
    const REAL one = 1.0;
    Cgemm("No transpose", "Conjugate transpose", k, m, n, COMPLEX(-one), a, lda, q, lda, COMPLEX(one), l, lda);
    //
    //     Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
    //
    REAL anorm = Clange("1", k, n, a, lda, rwork);
    REAL resid = Clange("1", k, m, l, lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max((INTEGER)1, n)).real()) / anorm) / eps;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute I - Q*Q'
    //
    Claset("Full", m, m, COMPLEX(zero), COMPLEX(one), l, lda);
    Cherk("Upper", "No transpose", m, n, -one, q, lda, one, l, lda);
    //
    //     Compute norm( I - Q*Q' ) / ( N * EPS ) .
    //
    resid = Clansy("1", "Upper", m, l, lda, rwork);
    //
    result[2 - 1] = (resid / (max((INTEGER)1, n)).real()) / eps;
    //
    //     End of Clqt02
    //
}
