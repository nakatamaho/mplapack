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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

void Cppt01(const char *uplo, INTEGER const n, COMPLEX *a, COMPLEX *afac, REAL *rwork, REAL &resid) {
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
    REAL anorm = Clanhp("1", uplo, n, a, rwork);
    const REAL one = 1.0;
    if (anorm <= zero) {
        resid = one / eps;
        return;
    }
    //
    //     Check the imaginary parts of the diagonal elements and return with
    //     an error code if any are nonzero.
    //
    INTEGER kc = 1;
    INTEGER k = 0;
    if (Mlsame(uplo, "U")) {
        for (k = 1; k <= n; k = k + 1) {
            if (afac[kc - 1].imag() != zero) {
                resid = one / eps;
                return;
            }
            kc += k + 1;
        }
    } else {
        for (k = 1; k <= n; k = k + 1) {
            if (afac[kc - 1].imag() != zero) {
                resid = one / eps;
                return;
            }
            kc += n - k + 1;
        }
    }
    //
    //     Compute the product U'*U, overwriting U.
    //
    REAL tr = 0.0;
    INTEGER i = 0;
    COMPLEX tc = 0.0;
    if (Mlsame(uplo, "U")) {
        kc = (n * (n - 1)) / 2 + 1;
        for (k = n; k >= 1; k = k - 1) {
            //
            //           Compute the (K,K) element of the result.
            //
            tr = Cdotc(k, afac[kc - 1], 1, afac[kc - 1], 1);
            afac[(kc + k - 1) - 1] = tr;
            //
            //           Compute the rest of column K.
            //
            if (k > 1) {
                Ctpmv("Upper", "Conjugate", "Non-unit", k - 1, afac, afac[kc - 1], 1);
                kc = kc - (k - 1);
            }
        }
        //
        //        Compute the difference  L*L' - A
        //
        kc = 1;
        for (k = 1; k <= n; k = k + 1) {
            for (i = 1; i <= k - 1; i = i + 1) {
                afac[(kc + i - 1) - 1] = afac[(kc + i - 1) - 1] - a[(kc + i - 1) - 1];
            }
            afac[(kc + k - 1) - 1] = afac[(kc + k - 1) - 1] - a[(kc + k - 1) - 1].real();
            kc += k;
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
                Chpr("Lower", n - k, one, afac[(kc + 1) - 1], 1, afac[(kc + n - k + 1) - 1]);
            }
            //
            //           Scale column K by the diagonal element.
            //
            tc = afac[kc - 1];
            Cscal(n - k + 1, tc, afac[kc - 1], 1);
            //
            kc = kc - (n - k + 2);
        }
        //
        //        Compute the difference  U'*U - A
        //
        kc = 1;
        for (k = 1; k <= n; k = k + 1) {
            afac[kc - 1] = afac[kc - 1] - a[kc - 1].real();
            for (i = k + 1; i <= n; i = i + 1) {
                afac[(kc + i - k) - 1] = afac[(kc + i - k) - 1] - a[(kc + i - k) - 1];
            }
            kc += n - k + 1;
        }
    }
    //
    //     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
    //
    resid = Clanhp("1", uplo, n, afac, rwork);
    //
    resid = ((resid / n.real()) / anorm) / eps;
    //
    //     End of Cppt01
    //
}
