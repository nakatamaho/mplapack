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

void Clsets(INTEGER const m, INTEGER const p, INTEGER const n, COMPLEX *a, COMPLEX *af, INTEGER const lda, COMPLEX *b, COMPLEX *bf, INTEGER const ldb, COMPLEX *c, COMPLEX *cf, COMPLEX *d, COMPLEX *df, COMPLEX *x, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result) {
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //
    //  ====================================================================
    //
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Copy the matrices A and B to the arrays AF and BF,
    //     and the vectors C and D to the arrays CF and DF,
    //
    Clacpy("Full", m, n, a, lda, af, lda);
    Clacpy("Full", p, n, b, ldb, bf, ldb);
    Ccopy(m, c, 1, cf, 1);
    Ccopy(p, d, 1, df, 1);
    //
    //     Solve LSE problem
    //
    INTEGER info = 0;
    Cgglse(m, n, p, af, lda, bf, ldb, cf, df, x, work, lwork, info);
    //
    //     Test the residual for the solution of LSE
    //
    //     Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS
    //
    Ccopy(m, c, 1, cf, 1);
    Ccopy(p, d, 1, df, 1);
    Cget02("No transpose", m, n, 1, a, lda, x, n, cf, m, rwork, result[1 - 1]);
    //
    //     Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS
    //
    Cget02("No transpose", p, n, 1, b, ldb, x, n, df, p, rwork, result[2 - 1]);
    //
    //     End of Clsets
    //
}
