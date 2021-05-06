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

REAL Rqrt17(const char *trans, INTEGER const iresid, INTEGER const m, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *x, INTEGER const ldx, REAL *b, INTEGER const ldb, REAL *c, REAL *work, INTEGER const lwork) {
    REAL return_value = 0.0;
    a([lda * star]);
    x([ldx * star]);
    b([ldb * star]);
    c([ldb * star]);
    work([lwork]);
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    return_value = zero;
    //
    INTEGER nrows = 0;
    INTEGER ncols = 0;
    if (Mlsame(trans, "N")) {
        nrows = m;
        ncols = n;
    } else if (Mlsame(trans, "T")) {
        nrows = n;
        ncols = m;
    } else {
        Mxerbla("Rqrt17", 1);
        return return_value;
    }
    //
    if (lwork < ncols * nrhs) {
        Mxerbla("Rqrt17", 13);
        return return_value;
    }
    //
    if (m <= 0 || n <= 0 || nrhs <= 0) {
        return return_value;
    }
    //
    REAL rwork[1];
    REAL norma = Rlange("One-norm", m, n, a, lda, rwork);
    REAL smlnum = Rlamch("Safe minimum") / Rlamch("Precision");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    INTEGER iscl = 0;
    //
    //     compute residual and scale it
    //
    Rlacpy("All", nrows, nrhs, b, ldb, c, ldb);
    Rgemm(trans, "No transpose", nrows, nrhs, ncols, -one, a, lda, x, ldx, one, c, ldb);
    REAL normrs = Rlange("Max", nrows, nrhs, c, ldb, rwork);
    INTEGER info = 0;
    if (normrs > smlnum) {
        iscl = 1;
        Rlascl("General", 0, 0, normrs, one, nrows, nrhs, c, ldb, info);
    }
    //
    //     compute R'*A
    //
    Rgemm("Transpose", trans, nrhs, ncols, nrows, one, c, ldb, a, lda, zero, work, nrhs);
    //
    //     compute and properly scale error
    //
    REAL err = Rlange("One-norm", nrhs, ncols, work, nrhs, rwork);
    if (norma != zero) {
        err = err / norma;
    }
    //
    if (iscl == 1) {
        err = err * normrs;
    }
    //
    REAL normb = 0.0;
    if (iresid == 1) {
        normb = Rlange("One-norm", nrows, nrhs, b, ldb, rwork);
        if (normb != zero) {
            err = err / normb;
        }
    } else {
        if (normrs != zero) {
            err = err / normrs;
        }
    }
    //
    return_value = err / (Rlamch("Epsilon") * castREAL(max({m, n, nrhs})));
    return return_value;
    //
    //     End of Rqrt17
    //
}
