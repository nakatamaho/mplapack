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

void Rppt01(const char *uplo, INTEGER const n, REAL *a, REAL *afac, REAL *rwork, REAL &resid) {
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
    if (n <= 0) {
        resid = zero;
        return;
    }
    //
    //     Exit with RESID = 1/EPS if ANORM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    REAL anorm = dlansp("1", uplo, n, a, rwork);
    const REAL one = 1.0;
    if (anorm <= zero) {
        resid = one / eps;
        return;
    }
    //
    //     Compute the product U'*U, overwriting U.
    //
    INTEGER kc = 0;
    INTEGER k = 0;
    REAL t = 0.0;
    if (Mlsame(uplo, "U")) {
        kc = (n * (n - 1)) / 2 + 1;
        for (k = n; k >= 1; k = k - 1) {
            //
            //           Compute the (K,K) element of the result.
            //
            t = Rdot(k, afac[kc - 1], 1, afac[kc - 1], 1);
            afac[(kc + k - 1) - 1] = t;
            //
            //           Compute the rest of column K.
            //
            if (k > 1) {
                Rtpmv("Upper", "Transpose", "Non-unit", k - 1, afac, afac[kc - 1], 1);
                kc = kc - (k - 1);
            }
        }
        //
        //     Compute the product L*L', overwriting L.
        //
    } else {
        kc = (n * (n + 1)) / 2;
        for (k = n; k >= 1; k = k - 1) {
            //
            //           Add a multiple of column K of the factor L to each of
            //           columns K+1 through N.
            //
            if (k < n) {
                Rspr("Lower", n - k, one, afac[(kc + 1) - 1], 1, afac[(kc + n - k + 1) - 1]);
            }
            //
            //           Scale column K by the diagonal element.
            //
            t = afac[kc - 1];
            Rscal(n - k + 1, t, afac[kc - 1], 1);
            //
            kc = kc - (n - k + 2);
        }
    }
    //
    //     Compute the difference  L*L' - A (or U'*U - A).
    //
    INTEGER npp = n * (n + 1) / 2;
    INTEGER i = 0;
    for (i = 1; i <= npp; i = i + 1) {
        afac[i - 1] = afac[i - 1] - a[i - 1];
    }
    //
    //     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
    //
    resid = dlansp("1", uplo, n, afac, rwork);
    //
    resid = ((resid / n.real()) / anorm) / eps;
    //
    //     End of Rppt01
    //
}
