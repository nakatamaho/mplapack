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

void Cppt03(const char *uplo, INTEGER const n, COMPLEX *a, COMPLEX *ainv, COMPLEX *work, INTEGER const ldwork, REAL *rwork, REAL &rcond, REAL &resid) {
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
    //     Quick exit if N = 0.
    //
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (n <= 0) {
        rcond = one;
        resid = zero;
        return;
    }
    //
    //     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    REAL anorm = Clanhp("1", uplo, n, a, rwork);
    REAL ainvnm = Clanhp("1", uplo, n, ainv, rwork);
    if (anorm <= zero || ainvnm <= zero) {
        rcond = zero;
        resid = one / eps;
        return;
    }
    rcond = (one / anorm) / ainvnm;
    //
    //     UPLO = 'U':
    //     Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
    //     expand it to a full matrix, then multiply by A one column at a
    //     time, moving the result one column to the left.
    //
    INTEGER jj = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (Mlsame(uplo, "U")) {
        //
        //        Copy AINV
        //
        jj = 1;
        for (j = 1; j <= n - 1; j = j + 1) {
            Ccopy(j, ainv[jj - 1], 1, &work[((j + 1) - 1) * ldwork], 1);
            for (i = 1; i <= j - 1; i = i + 1) {
                work[(j - 1) + ((i + 1) - 1) * ldwork] = conj(ainv[(jj + i - 1) - 1]);
            }
            jj += j;
        }
        jj = ((n - 1) * n) / 2 + 1;
        for (i = 1; i <= n - 1; i = i + 1) {
            work[(n - 1) + ((i + 1) - 1) * ldwork] = conj(ainv[(jj + i - 1) - 1]);
        }
        //
        //        Multiply by A
        //
        for (j = 1; j <= n - 1; j = j + 1) {
            Chpmv("Upper", n, -cone, a, &work[((j + 1) - 1) * ldwork], 1, czero, &work[(j - 1) * ldwork], 1);
        }
        Chpmv("Upper", n, -cone, a, ainv[jj - 1], 1, czero, &work[(n - 1) * ldwork], 1);
        //
        //     UPLO = 'L':
        //     Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
        //     and multiply by A, moving each column to the right.
        //
    } else {
        //
        //        Copy AINV
        //
        for (i = 1; i <= n - 1; i = i + 1) {
            work[(i - 1) * ldwork] = conj(ainv[(i + 1) - 1]);
        }
        jj = n + 1;
        for (j = 2; j <= n; j = j + 1) {
            Ccopy(n - j + 1, ainv[jj - 1], 1, &work[(j - 1) + ((j - 1) - 1) * ldwork], 1);
            for (i = 1; i <= n - j; i = i + 1) {
                work[(j - 1) + ((j + i - 1) - 1) * ldwork] = conj(ainv[(jj + i) - 1]);
            }
            jj += n - j + 1;
        }
        //
        //        Multiply by A
        //
        for (j = n; j >= 2; j = j - 1) {
            Chpmv("Lower", n, -cone, a, &work[((j - 1) - 1) * ldwork], 1, czero, &work[(j - 1) * ldwork], 1);
        }
        Chpmv("Lower", n, -cone, a, ainv[1 - 1], 1, czero, &work[(1 - 1)], 1);
        //
    }
    //
    //     Add the identity matrix to WORK .
    //
    for (i = 1; i <= n; i = i + 1) {
        work[(i - 1) + (i - 1) * ldwork] += cone;
    }
    //
    //     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
    //
    resid = Clange("1", n, n, work, ldwork, rwork);
    //
    resid = ((resid * rcond) / eps) / n.real();
    //
    //     End of Cppt03
    //
}
