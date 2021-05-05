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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Rgrqts(INTEGER const m, INTEGER const p, INTEGER const n, REAL *a, REAL *af, REAL *q, REAL *r, INTEGER const lda, REAL *taua, REAL *b, REAL *bf, REAL *z, REAL *t, REAL *bwk, INTEGER const ldb, REAL *taub, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
    a([lda * star]);
    af([lda * star]);
    q([lda * star]);
    r([lda * star]);
    b([ldb * star]);
    bf([ldb * star]);
    z([ldb * star]);
    t([ldb * star]);
    bwk([ldb * star]);
    work([lwork]);
    result([4]);
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
    REAL ulp = Rlamch("Precision");
    REAL unfl = Rlamch("Safe minimum");
    //
    //     Copy the matrix A to the array AF.
    //
    Rlacpy("Full", m, n, a, lda, af, lda);
    Rlacpy("Full", p, n, b, ldb, bf, ldb);
    //
    REAL anorm = max({Rlange("1", m, n, a, lda, rwork), unfl});
    REAL bnorm = max({Rlange("1", p, n, b, ldb, rwork), unfl});
    //
    //     Factorize the matrices A and B in the arrays AF and BF.
    //
    INTEGER info = 0;
    Rggrqf(m, p, n, af, lda, taua, bf, ldb, taub, work, lwork, info);
    //
    //     Generate the N-by-N matrix Q
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
    Rorgrq(n, n, min(m, n), q, lda, taua, work, lwork, info);
    //
    //     Generate the P-by-P matrix Z
    //
    Rlaset("Full", p, p, rogue, rogue, z, ldb);
    if (p > 1) {
        Rlacpy("Lower", p - 1, n, bf[(2 - 1)], ldb, &z[(2 - 1)], ldb);
    }
    Rorgqr(p, p, min(p, n), z, ldb, taub, work, lwork, info);
    //
    //     Copy R
    //
    const REAL zero = 0.0;
    Rlaset("Full", m, n, zero, zero, r, lda);
    if (m <= n) {
        Rlacpy("Upper", m, m, af[((n - m + 1) - 1) * ldaf], lda, r[((n - m + 1) - 1) * ldr], lda);
    } else {
        Rlacpy("Full", m - n, n, af, lda, r, lda);
        Rlacpy("Upper", n, n, af[((m - n + 1) - 1)], lda, r[((m - n + 1) - 1)], lda);
    }
    //
    //     Copy T
    //
    Rlaset("Full", p, n, zero, zero, t, ldb);
    Rlacpy("Upper", p, n, bf, ldb, t, ldb);
    //
    //     Compute R - A*Q'
    //
    const REAL one = 1.0;
    Rgemm("No transpose", "Transpose", m, n, n, -one, a, lda, q, lda, one, r, lda);
    //
    //     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
    //
    REAL resid = Rlange("1", m, n, r, lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max({(INTEGER)1, m, n})).real()) / anorm) / ulp;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute T*Q - Z'*B
    //
    Rgemm("Transpose", "No transpose", p, n, p, one, z, ldb, b, ldb, zero, bwk, ldb);
    Rgemm("No transpose", "No transpose", p, n, n, one, t, ldb, q, lda, -one, bwk, ldb);
    //
    //     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
    //
    resid = Rlange("1", p, n, bwk, ldb, rwork);
    if (bnorm > zero) {
        result[2 - 1] = ((resid / (max({(INTEGER)1, p, m})).real()) / bnorm) / ulp;
    } else {
        result[2 - 1] = zero;
    }
    //
    //     Compute I - Q*Q'
    //
    Rlaset("Full", n, n, zero, one, r, lda);
    Rsyrk("Upper", "No Transpose", n, n, -one, q, lda, one, r, lda);
    //
    //     Compute norm( I - Q'*Q ) / ( N * ULP ) .
    //
    resid = Rlansy("1", "Upper", n, r, lda, rwork);
    result[3 - 1] = (resid / (max((INTEGER)1, n)).real()) / ulp;
    //
    //     Compute I - Z'*Z
    //
    Rlaset("Full", p, p, zero, one, t, ldb);
    Rsyrk("Upper", "Transpose", p, p, -one, z, ldb, one, t, ldb);
    //
    //     Compute norm( I - Z'*Z ) / ( P*ULP ) .
    //
    resid = Rlansy("1", "Upper", p, t, ldb, rwork);
    result[4 - 1] = (resid / (max((INTEGER)1, p)).real()) / ulp;
    //
    //     End of Rgrqts
    //
}
