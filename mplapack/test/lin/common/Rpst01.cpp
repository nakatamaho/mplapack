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

void Rpst01(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *afac, INTEGER const ldafac, REAL *perm, INTEGER const ldperm, INTEGER *piv, REAL *rwork, REAL &resid, INTEGER const rank) {
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
    //     Quick exit if N = 0.
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
    REAL anorm = Rlansy("1", uplo, n, a, lda, rwork);
    const REAL one = 1.0;
    if (anorm <= zero) {
        resid = one / eps;
        return;
    }
    //
    //     Compute the product U'*U, overwriting U.
    //
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER k = 0;
    REAL t = 0.0;
    if (Mlsame(uplo, "U")) {
        //
        if (rank < n) {
            for (j = rank + 1; j <= n; j = j + 1) {
                for (i = rank + 1; i <= j; i = i + 1) {
                    afac[(i - 1) + (j - 1) * ldafac] = zero;
                }
            }
        }
        //
        for (k = n; k >= 1; k = k - 1) {
            //
            //           Compute the (K,K) element of the result.
            //
            t = Rdot(k, &afac[(k - 1) * ldafac], 1, &afac[(k - 1) * ldafac], 1);
            afac[(k - 1) + (k - 1) * ldafac] = t;
            //
            //           Compute the rest of column K.
            //
            Rtrmv("Upper", "Transpose", "Non-unit", k - 1, afac, ldafac, &afac[(k - 1) * ldafac], 1);
            //
        }
        //
        //     Compute the product L*L', overwriting L.
        //
    } else {
        //
        if (rank < n) {
            for (j = rank + 1; j <= n; j = j + 1) {
                for (i = j; i <= n; i = i + 1) {
                    afac[(i - 1) + (j - 1) * ldafac] = zero;
                }
            }
        }
        //
        for (k = n; k >= 1; k = k - 1) {
            //           Add a multiple of column K of the factor L to each of
            //           columns K+1 through N.
            //
            if (k + 1 <= n) {
                Rsyr("Lower", n - k, one, &afac[((k + 1) - 1) + (k - 1) * ldafac], 1, &afac[((k + 1) - 1) + ((k + 1) - 1) * ldafac], ldafac);
            }
            //
            //           Scale column K by the diagonal element.
            //
            t = afac[(k - 1) + (k - 1) * ldafac];
            Rscal(n - k + 1, t, &afac[(k - 1) + (k - 1) * ldafac], 1);
        }
        //
    }
    //
    //        Form P*L*L'*P' or P*U'*U*P'
    //
    if (Mlsame(uplo, "U")) {
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                if (piv[i - 1] <= piv[j - 1]) {
                    if (i <= j) {
                        perm[(piv[i - 1] - 1) + (piv[j - 1] - 1) * ldperm] = afac[(i - 1) + (j - 1) * ldafac];
                    } else {
                        perm[(piv[i - 1] - 1) + (piv[j - 1] - 1) * ldperm] = afac[(j - 1) + (i - 1) * ldafac];
                    }
                }
            }
        }
        //
    } else {
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                if (piv[i - 1] >= piv[j - 1]) {
                    if (i >= j) {
                        perm[(piv[i - 1] - 1) + (piv[j - 1] - 1) * ldperm] = afac[(i - 1) + (j - 1) * ldafac];
                    } else {
                        perm[(piv[i - 1] - 1) + (piv[j - 1] - 1) * ldperm] = afac[(j - 1) + (i - 1) * ldafac];
                    }
                }
            }
        }
        //
    }
    //
    //     Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).
    //
    if (Mlsame(uplo, "U")) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                perm[(i - 1) + (j - 1) * ldperm] = perm[(i - 1) + (j - 1) * ldperm] - a[(i - 1) + (j - 1) * lda];
            }
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            for (i = j; i <= n; i = i + 1) {
                perm[(i - 1) + (j - 1) * ldperm] = perm[(i - 1) + (j - 1) * ldperm] - a[(i - 1) + (j - 1) * lda];
            }
        }
    }
    //
    //     Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
    //     ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).
    //
    resid = Rlansy("1", uplo, n, perm, ldafac, rwork);
    //
    resid = ((resid / castREAL(n)) / anorm) / eps;
    //
    //     End of Rpst01
    //
}
