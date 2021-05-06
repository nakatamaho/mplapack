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

void Rppt03(const char *uplo, INTEGER const n, REAL *a, REAL *ainv, REAL *work, INTEGER const ldwork, REAL *rwork, REAL &rcond, REAL &resid) {
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
    REAL anorm = Rlansp("1", uplo, n, a, rwork);
    REAL ainvnm = Rlansp("1", uplo, n, ainv, rwork);
    if (anorm <= zero || ainvnm == zero) {
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
    if (Mlsame(uplo, "U")) {
        //
        //        Copy AINV
        //
        jj = 1;
        for (j = 1; j <= n - 1; j = j + 1) {
            Rcopy(j, &ainv[jj - 1], 1, &work[((j + 1) - 1) * ldwork], 1);
            Rcopy(j - 1, &ainv[jj - 1], 1, &work[(j - 1) + (2 - 1) * ldwork], ldwork);
            jj += j;
        }
        jj = ((n - 1) * n) / 2 + 1;
        Rcopy(n - 1, &ainv[jj - 1], 1, &work[(n - 1) + (2 - 1) * ldwork], ldwork);
        //
        //        Multiply by A
        //
        for (j = 1; j <= n - 1; j = j + 1) {
            Rspmv("Upper", n, -one, a, &work[((j + 1) - 1) * ldwork], 1, zero, &work[(j - 1) * ldwork], 1);
        }
        Rspmv("Upper", n, -one, a, &ainv[jj - 1], 1, zero, &work[(n - 1) * ldwork], 1);
        //
        //     UPLO = 'L':
        //     Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
        //     and multiply by A, moving each column to the right.
        //
    } else {
        //
        //        Copy AINV
        //
        Rcopy(n - 1, &ainv[2 - 1], 1, &work[(1 - 1)], ldwork);
        jj = n + 1;
        for (j = 2; j <= n; j = j + 1) {
            Rcopy(n - j + 1, &ainv[jj - 1], 1, &work[(j - 1) + ((j - 1) - 1) * ldwork], 1);
            Rcopy(n - j, &ainv[(jj + 1) - 1], 1, &work[(j - 1) + (j - 1) * ldwork], ldwork);
            jj += n - j + 1;
        }
        //
        //        Multiply by A
        //
        for (j = n; j >= 2; j = j - 1) {
            Rspmv("Lower", n, -one, a, &work[((j - 1) - 1) * ldwork], 1, zero, &work[(j - 1) * ldwork], 1);
        }
        Rspmv("Lower", n, -one, a, &ainv[1 - 1], 1, zero, &work[(1 - 1)], 1);
        //
    }
    //
    //     Add the identity matrix to WORK .
    //
    INTEGER i = 0;
    for (i = 1; i <= n; i = i + 1) {
        work[(i - 1) + (i - 1) * ldwork] += one;
    }
    //
    //     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
    //
    resid = Rlange("1", n, n, work, ldwork, rwork);
    //
    resid = ((resid * rcond) / eps) / castREAL(n);
    //
    //     End of Rppt03
    //
}
