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

void Rgtt01(INTEGER const n, REAL *dl, REAL *d, REAL *du, REAL *dlf, REAL *df, REAL *duf, REAL *du2, INTEGER *ipiv, REAL *work, INTEGER const ldwork, REAL *rwork, REAL &resid) {
    work([ldwork * star]);
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (n <= 0) {
        resid = zero;
        return;
    }
    //
    REAL eps = Rlamch("Epsilon");
    //
    //     Copy the matrix U to WORK.
    //
    INTEGER j = 0;
    INTEGER i = 0;
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= n; i = i + 1) {
            work[(i - 1) + (j - 1) * ldwork] = zero;
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        if (i == 1) {
            work[(i - 1) + (i - 1) * ldwork] = df[i - 1];
            if (n >= 2) {
                work[(i - 1) + ((i + 1) - 1) * ldwork] = duf[i - 1];
            }
            if (n >= 3) {
                work[(i - 1) + ((i + 2) - 1) * ldwork] = du2[i - 1];
            }
        } else if (i == n) {
            work[(i - 1) + (i - 1) * ldwork] = df[i - 1];
        } else {
            work[(i - 1) + (i - 1) * ldwork] = df[i - 1];
            work[(i - 1) + ((i + 1) - 1) * ldwork] = duf[i - 1];
            if (i < n - 1) {
                work[(i - 1) + ((i + 2) - 1) * ldwork] = du2[i - 1];
            }
        }
    }
    //
    //     Multiply on the left by L.
    //
    INTEGER lastj = n;
    REAL li = 0.0;
    INTEGER ip = 0;
    for (i = n - 1; i >= 1; i = i - 1) {
        li = dlf[i - 1];
        Raxpy(lastj - i + 1, li, &work[(i - 1) + (i - 1) * ldwork], ldwork, &work[((i + 1) - 1) + (i - 1) * ldwork], ldwork);
        ip = ipiv[i - 1];
        if (ip == i) {
            lastj = min(i + 2, n);
        } else {
            Rswap(lastj - i + 1, &work[(i - 1) + (i - 1) * ldwork], ldwork, &work[((i + 1) - 1) + (i - 1) * ldwork], ldwork);
        }
    }
    //
    //     Subtract the matrix A.
    //
    work[(1 - 1)] = work[(1 - 1)] - d[1 - 1];
    if (n > 1) {
        work[(2 - 1) * ldwork] = work[(2 - 1) * ldwork] - du[1 - 1];
        work[(n - 1) + ((n - 1) - 1) * ldwork] = work[(n - 1) + ((n - 1) - 1) * ldwork] - dl[(n - 1) - 1];
        work[(n - 1) + (n - 1) * ldwork] = work[(n - 1) + (n - 1) * ldwork] - d[n - 1];
        for (i = 2; i <= n - 1; i = i + 1) {
            work[(i - 1) + ((i - 1) - 1) * ldwork] = work[(i - 1) + ((i - 1) - 1) * ldwork] - dl[(i - 1) - 1];
            work[(i - 1) + (i - 1) * ldwork] = work[(i - 1) + (i - 1) * ldwork] - d[i - 1];
            work[(i - 1) + ((i + 1) - 1) * ldwork] = work[(i - 1) + ((i + 1) - 1) * ldwork] - du[i - 1];
        }
    }
    //
    //     Compute the 1-norm of the tridiagonal matrix A.
    //
    REAL anorm = Rlangt("1", n, dl, d, du);
    //
    //     Compute the 1-norm of WORK, which is only guaranteed to be
    //     upper Hessenberg.
    //
    resid = Rlanhs("1", n, work, ldwork, rwork);
    //
    //     Compute norm(L*U - A) / (norm(A) * EPS)
    //
    const REAL one = 1.0;
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        resid = (resid / anorm) / eps;
    }
    //
    //     End of Rgtt01
    //
}
