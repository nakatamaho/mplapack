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

// Derived from LAPACK routine ZQRT17.
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

REAL Cqrt17(fem::str_cref trans, INTEGER const iresid, INTEGER const m, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, COMPLEX *c, COMPLEX *work, INTEGER const lwork) {
    REAL return_value = 0.0;
    //
    const REAL zero = 0.0;
    return_value = zero;
    //
    INTEGER nrows = 0;
    INTEGER ncols = 0;
    if (Mlsame(trans.elems(), "N")) {
        nrows = m;
        ncols = n;
    } else if (Mlsame(trans.elems(), "C")) {
        nrows = n;
        ncols = m;
    } else {
        Mxerbla("Cqrt17", 1);
        return return_value;
    }
    //
    if (lwork < ncols * nrhs) {
        Mxerbla("Cqrt17", 13);
        return return_value;
    }
    //
    if (m <= 0 || n <= 0 || nrhs <= 0) {
        return return_value;
    }
    //
    REAL rwork[1];
    REAL norma = Clange("One-norm", m, n, a, lda, rwork);
    REAL smlnum = Rlamch("Safe minimum") / Rlamch("Precision");
    INTEGER iscl = 0;
    //
    // compute residual and scale it
    //
    Clacpy("All", nrows, nrhs, b, ldb, c, ldb);
    const REAL one = 1.0;
    Cgemm(trans.elems(), "No transpose", nrows, nrhs, ncols, COMPLEX(-one), a, lda, x, ldx, COMPLEX(one), c, ldb);
    REAL normrs = Clange("Max", nrows, nrhs, c, ldb, rwork);
    INTEGER info = 0;
    if (normrs > smlnum) {
        iscl = 1;
        Clascl("General", 0, 0, normrs, one, nrows, nrhs, c, ldb, info);
    }
    //
    // compute R**H * op(A)
    //
    Cgemm("Conjugate transpose", trans.elems(), nrhs, ncols, nrows, COMPLEX(one), c, ldb, a, lda, COMPLEX(zero), work, nrhs);
    //
    // compute and properly scale error
    //
    REAL err = Clange("One-norm", nrhs, ncols, work, nrhs, rwork);
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
        normb = Clange("One-norm", nrows, nrhs, b, ldb, rwork);
        if (normb != zero) {
            err = err / normb;
        }
    } else {
        if (normrs != zero) {
            err = err / normrs;
        }
    }
    //
    return_value = err / (Rlamch("Epsilon") * castREAL(max(m, n, nrhs)));
    return return_value;
    //
    // End of Cqrt17
    //
}
