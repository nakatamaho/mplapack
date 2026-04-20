/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DTRT02.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Rtrt02(fem::str_cref uplo, fem::str_cref trans, fem::str_cref diag, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *x, INTEGER const ldx, REAL *b, INTEGER const ldb, REAL *work, REAL &resid) {
    //
    // Quick exit if N = 0 or NRHS = 0
    //
    const REAL zero = 0.0;
    if (n <= 0 || nrhs <= 0) {
        resid = zero;
        return;
    }
    //
    // Compute the 1-norm of op(A).
    //
    REAL anorm = 0.0;
    if (Mlsame(trans.elems(), "N")) {
        anorm = Rlantr("1", uplo.elems(), diag.elems(), n, n, a, lda, work);
    } else {
        anorm = Rlantr("I", uplo.elems(), diag.elems(), n, n, a, lda, work);
    }
    //
    // Exit with RESID = 1/EPS if ANORM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    const REAL one = 1.0;
    if (anorm <= zero) {
        resid = one / eps;
        return;
    }
    //
    // Compute the maximum over the number of right hand sides of
    // norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS )
    //
    resid = zero;
    INTEGER j = 0;
    REAL bnorm = 0.0;
    REAL xnorm = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        Rcopy(n, &x[(j - 1) * ldx], 1, work, 1);
        Rtrmv(uplo.elems(), trans.elems(), diag.elems(), n, a, lda, work, 1);
        Raxpy(n, -one, &b[(j - 1) * ldb], 1, work, 1);
        bnorm = Rasum(n, work, 1);
        xnorm = Rasum(n, &x[(j - 1) * ldx], 1);
        if (xnorm <= zero) {
            resid = one / eps;
        } else {
            resid = max(resid, ((bnorm / anorm) / xnorm) / eps);
        }
    }
    //
    // End of Rtrt02
    //
}
