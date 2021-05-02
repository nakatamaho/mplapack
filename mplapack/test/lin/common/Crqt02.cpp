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

void Crqt02(INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, COMPLEX *af, COMPLEX *q, COMPLEX *r, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    //     Copy the last k rows of the factorization to the array Q
    //
    const COMPLEX rogue = COMPLEX(-1.0e+10, -1.0e+10);
    Claset("Full", m, n, rogue, rogue, q, lda);
    if (k < n) {
        Clacpy("Full", k, n - k, af[((m - k + 1) - 1)], lda, &q[((m - k + 1) - 1)], lda);
    }
    if (k > 1) {
        Clacpy("Lower", k - 1, k - 1, af[((m - k + 2) - 1) + ((n - k + 1) - 1) * ldaf], lda, &q[((m - k + 2) - 1) + ((n - k + 1) - 1) * ldq], lda);
    }
    //
    //     Generate the last n rows of the matrix Q
    //
    cmn.srnamt = "Cungrq";
    INTEGER info = 0;
    Cungrq(m, n, k, q, lda, &tau[(m - k + 1) - 1], work, lwork, info);
    //
    //     Copy R(m-k+1:m,n-m+1:n)
    //
    Claset("Full", k, m, COMPLEX(zero), COMPLEX(zero), r[((m - k + 1) - 1) + ((n - m + 1) - 1) * ldr], lda);
    Clacpy("Upper", k, k, af[((m - k + 1) - 1) + ((n - k + 1) - 1) * ldaf], lda, r[((m - k + 1) - 1) + ((n - k + 1) - 1) * ldr], lda);
    //
    //     Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'
    //
    const REAL one = 1.0;
    Cgemm("No transpose", "Conjugate transpose", k, m, n, COMPLEX(-one), &a[((m - k + 1) - 1)], lda, q, lda, COMPLEX(one), r[((m - k + 1) - 1) + ((n - m + 1) - 1) * ldr], lda);
    //
    //     Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .
    //
    REAL anorm = Clange("1", k, n, &a[((m - k + 1) - 1)], lda, rwork);
    REAL resid = Clange("1", k, m, r[((m - k + 1) - 1) + ((n - m + 1) - 1) * ldr], lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max((INTEGER)1, n)).real()) / anorm) / eps;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute I - Q*Q'
    //
    Claset("Full", m, m, COMPLEX(zero), COMPLEX(one), r, lda);
    Cherk("Upper", "No transpose", m, n, -one, q, lda, one, r, lda);
    //
    //     Compute norm( I - Q*Q' ) / ( N * EPS ) .
    //
    resid = Clansy("1", "Upper", m, r, lda, rwork);
    //
    result[2 - 1] = (resid / (max((INTEGER)1, n)).real()) / eps;
    //
    //     End of Crqt02
    //
}
