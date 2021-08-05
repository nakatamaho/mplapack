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

void Rhst01(INTEGER const n, INTEGER const /* ilo */, INTEGER const /* ihi */, REAL *a, INTEGER const lda, REAL *h, INTEGER const ldh, REAL *q, INTEGER const ldq, REAL *work, INTEGER const lwork, REAL *result) {
    //
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (n <= 0) {
        result[1 - 1] = zero;
        result[2 - 1] = zero;
        return;
    }
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL eps = Rlamch("Precision");
    const REAL one = 1.0;
    REAL ovfl = one / unfl;
    REAL smlnum = unfl * castREAL(n) / eps;
    //
    //     Test 1:  Compute norm( A - Q*H*Q' ) / ( norm(A) * N * EPS )
    //
    //     Copy A to WORK
    //
    INTEGER ldwork = max((INTEGER)1, n);
    Rlacpy(" ", n, n, a, lda, work, ldwork);
    //
    //     Compute Q*H
    //
    Rgemm("No transpose", "No transpose", n, n, n, one, q, ldq, h, ldh, zero, &work[(ldwork * n + 1) - 1], ldwork);
    //
    //     Compute A - Q*H*Q'
    //
    Rgemm("No transpose", "Transpose", n, n, n, -one, &work[(ldwork * n + 1) - 1], ldwork, q, ldq, one, work, ldwork);
    //
    REAL anorm = max({Rlange("1", n, n, a, lda, &work[(ldwork * n + 1) - 1]), unfl});
    REAL wnorm = Rlange("1", n, n, work, ldwork, &work[(ldwork * n + 1) - 1]);
    //
    //     Note that RESULT(1) cannot overflow and is bounded by 1/(N*EPS)
    //
    result[1 - 1] = min(wnorm, anorm) / max(smlnum, REAL(anorm * eps)) / castREAL(n);
    //
    //     Test 2:  Compute norm( I - Q'*Q ) / ( N * EPS )
    //
    Rort01("Columns", n, n, q, ldq, work, lwork, result[2 - 1]);
    //
    //     End of Rhst01
    //
}
