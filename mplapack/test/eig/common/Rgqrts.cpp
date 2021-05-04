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

void Rgqrts(INTEGER const n, INTEGER const m, INTEGER const p, REAL *a, REAL *af, REAL *q, REAL *r, INTEGER const lda, REAL *taua, REAL *b, REAL *bf, REAL *z, REAL *t, REAL *bwk, INTEGER const ldb, REAL *taub, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    Rlacpy("Full", n, m, a, lda, af, lda);
    Rlacpy("Full", n, p, b, ldb, bf, ldb);
    //
    REAL anorm = max({Rlange("1", n, m, a, lda, rwork), unfl});
    REAL bnorm = max({Rlange("1", n, p, b, ldb, rwork), unfl});
    //
    //     Factorize the matrices A and B in the arrays AF and BF.
    //
    INTEGER info = 0;
    Rggqrf(n, m, p, af, lda, taua, bf, ldb, taub, work, lwork, info);
    //
    //     Generate the N-by-N matrix Q
    //
    const REAL rogue = -1.0e+10;
    Rlaset("Full", n, n, rogue, rogue, q, lda);
    Rlacpy("Lower", n - 1, m, af[(2 - 1)], lda, &q[(2 - 1)], lda);
    Rorgqr(n, n, min(n, m), q, lda, taua, work, lwork, info);
    //
    //     Generate the P-by-P matrix Z
    //
    Rlaset("Full", p, p, rogue, rogue, z, ldb);
    if (n <= p) {
        if (n > 0 && n < p) {
            Rlacpy("Full", n, p - n, bf, ldb, &z[((p - n + 1) - 1)], ldb);
        }
        if (n > 1) {
            Rlacpy("Lower", n - 1, n - 1, bf[(2 - 1) + ((p - n + 1) - 1) * ldbf], ldb, &z[((p - n + 2) - 1) + ((p - n + 1) - 1) * ldz], ldb);
        }
    } else {
        if (p > 1) {
            Rlacpy("Lower", p - 1, p - 1, bf[((n - p + 2) - 1)], ldb, &z[(2 - 1)], ldb);
        }
    }
    Rorgrq(p, p, min(n, p), z, ldb, taub, work, lwork, info);
    //
    //     Copy R
    //
    const REAL zero = 0.0;
    Rlaset("Full", n, m, zero, zero, r, lda);
    Rlacpy("Upper", n, m, af, lda, r, lda);
    //
    //     Copy T
    //
    Rlaset("Full", n, p, zero, zero, t, ldb);
    if (n <= p) {
        Rlacpy("Upper", n, n, bf[((p - n + 1) - 1) * ldbf], ldb, &t[((p - n + 1) - 1) * ldt], ldb);
    } else {
        Rlacpy("Full", n - p, p, bf, ldb, t, ldb);
        Rlacpy("Upper", p, p, bf[((n - p + 1) - 1)], ldb, &t[((n - p + 1) - 1)], ldb);
    }
    //
    //     Compute R - Q'*A
    //
    const REAL one = 1.0;
    Rgemm("Transpose", "No transpose", n, m, n, -one, q, lda, a, lda, one, r, lda);
    //
    //     Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .
    //
    REAL resid = Rlange("1", n, m, r, lda, rwork);
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max({(INTEGER)1, m, n})).real()) / anorm) / ulp;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute T*Z - Q'*B
    //
    Rgemm("No Transpose", "No transpose", n, p, p, one, t, ldb, z, ldb, zero, bwk, ldb);
    Rgemm("Transpose", "No transpose", n, p, n, -one, q, lda, b, ldb, one, bwk, ldb);
    //
    //     Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
    //
    resid = Rlange("1", n, p, bwk, ldb, rwork);
    if (bnorm > zero) {
        result[2 - 1] = ((resid / (max({(INTEGER)1, p, n})).real()) / bnorm) / ulp;
    } else {
        result[2 - 1] = zero;
    }
    //
    //     Compute I - Q'*Q
    //
    Rlaset("Full", n, n, zero, one, r, lda);
    Rsyrk("Upper", "Transpose", n, n, -one, q, lda, one, r, lda);
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
    //     End of Rgqrts
    //
}
