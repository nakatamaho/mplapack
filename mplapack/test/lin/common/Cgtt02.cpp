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
#include <mplapack.h>

void Cgtt02(const char *trans, INTEGER const n, INTEGER const nrhs, COMPLEX *dl, COMPLEX *d, COMPLEX *du, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, REAL &resid) {
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
    //     Quick exit if N = 0 or NRHS = 0
    //
    const REAL zero = 0.0;
    resid = zero;
    if (n <= 0 || nrhs == 0) {
        return;
    }
    //
    //     Compute the maximum over the number of right hand sides of
    //        norm(B - op(A)*X) / ( norm(A) * norm(X) * EPS ).
    //
    REAL anorm = 0.0;
    if (Mlsame(trans, "N")) {
        anorm = zlangt("1", n, dl, d, du);
    } else {
        anorm = zlangt("I", n, dl, d, du);
    }
    //
    //     Exit with RESID = 1/EPS if ANORM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    const REAL one = 1.0;
    if (anorm <= zero) {
        resid = one / eps;
        return;
    }
    //
    //     Compute B - op(A)*X.
    //
    zlagtm(trans, n, nrhs, -one, dl, d, du, x, ldx, one, b, ldb);
    //
    INTEGER j = 0;
    REAL bnorm = 0.0;
    REAL xnorm = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        bnorm = RCasum(n, &b[(j - 1) * ldb], 1);
        xnorm = RCasum(n, &x[(j - 1) * ldx], 1);
        if (xnorm <= zero) {
            resid = one / eps;
        } else {
            resid = max(resid, ((bnorm / anorm) / xnorm) / eps);
        }
    }
    //
    //     End of Cgtt02
    //
}
