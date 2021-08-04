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

void Ctbt03(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const kd, INTEGER const nrhs, COMPLEX *ab, INTEGER const ldab, REAL const scale, REAL *cnorm, REAL const tscal, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, COMPLEX *work, REAL &resid) {
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
    //     Quick exit if N = 0
    //
    const REAL zero = 0.0;
    if (n <= 0 || nrhs <= 0) {
        resid = zero;
        return;
    }
    REAL eps = Rlamch("Epsilon");
    REAL smlnum = Rlamch("Safe minimum");
    //
    //     Compute the norm of the triangular matrix A using the column
    //     norms already computed by Clatbs.
    //
    REAL tnorm = zero;
    INTEGER j = 0;
    if (Mlsame(diag, "N")) {
        if (Mlsame(uplo, "U")) {
            for (j = 1; j <= n; j = j + 1) {
                tnorm = max(tnorm, REAL(tscal * abs(ab[((kd + 1) - 1) + (j - 1) * ldab]) + cnorm[j - 1]));
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                tnorm = max(tnorm, REAL(tscal * abs(ab[(j - 1) * ldab]) + cnorm[j - 1]));
            }
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            tnorm = max(tnorm, REAL(tscal + cnorm[j - 1]));
        }
    }
    //
    //     Compute the maximum over the number of right hand sides of
    //        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
    //
    resid = zero;
    INTEGER ix = 0;
    const REAL one = 1.0;
    REAL xnorm = 0.0;
    REAL xscal = 0.0;
    REAL err = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        Ccopy(n, &x[(j - 1) * ldx], 1, work, 1);
        ix = iCamax(n, work, 1);
        xnorm = max(one, abs(x[(ix - 1) + (j - 1) * ldx]));
        xscal = (one / xnorm) / castREAL(kd + 1);
        CRscal(n, xscal, work, 1);
        Ctbmv(uplo, trans, diag, n, kd, ab, ldab, work, 1);
        Caxpy(n, COMPLEX(-scale * xscal), &b[(j - 1) * ldb], 1, work, 1);
        ix = iCamax(n, work, 1);
        err = tscal * abs(work[ix - 1]);
        ix = iCamax(n, &x[(j - 1) * ldx], 1);
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
    //     End of Ctbt03
    //
}
