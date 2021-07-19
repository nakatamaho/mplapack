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

void Rtrt03(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL const scale, REAL *cnorm, REAL const tscal, REAL *x, INTEGER const ldx, REAL *b, INTEGER const ldb, REAL *work, REAL &resid) {
    //
    //     Quick exit if N = 0
    //
    const REAL zero = 0.0;
    if (n <= 0 || nrhs <= 0) {
        resid = zero;
        return;
    }
    REAL eps = Rlamch("Epsilon");
    REAL smlnum = Rlamch("Safe minimum");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    //     Compute the norm of the triangular matrix A using the column
    //     norms already computed by Rlatrs.
    //
    REAL tnorm = zero;
    INTEGER j = 0;
    if (Mlsame(diag, "N")) {
        for (j = 1; j <= n; j = j + 1) {
            tnorm = max(tnorm, tscal * abs(a[(j - 1) + (j - 1) * lda]) + cnorm[j - 1]);
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            tnorm = max(tnorm, tscal + cnorm[j - 1]);
        }
    }
    //
    //     Compute the maximum over the number of right hand sides of
    //        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
    //
    resid = zero;
    INTEGER ix = 0;
    REAL xnorm = 0.0;
    REAL xscal = 0.0;
    REAL err = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        Rcopy(n, &x[(j - 1) * ldx], 1, work, 1);
        ix = iRamax(n, work, 1);
        xnorm = max(one, abs(x[(ix - 1) + (j - 1) * ldx]));
        xscal = (one / xnorm) / castREAL(n);
        Rscal(n, xscal, work, 1);
        Rtrmv(uplo, trans, diag, n, a, lda, work, 1);
        Raxpy(n, -scale * xscal, &b[(j - 1) * ldb], 1, work, 1);
        ix = iRamax(n, work, 1);
        err = tscal * abs(work[ix - 1]);
        ix = iRamax(n, &x[(j - 1) * ldx], 1);
        xnorm = abs(x[(ix - 1) + (j - 1) * ldx]);
        if (err * smlnum <= xnorm) {
            if (xnorm > zero) {
                err = err / xnorm;
            }
        } else {
            if (err > zero) {
                err = one / eps;
            }
        }
        if (err * smlnum <= tnorm) {
            if (tnorm > zero) {
                err = err / tnorm;
            }
        } else {
            if (err > zero) {
                err = one / eps;
            }
        }
        resid = max(resid, err);
    }
    //
    //     End of Rtrt03
    //
}
