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

void Cspt03(const char *uplo, INTEGER const n, COMPLEX *a, COMPLEX *ainv, COMPLEX *work, INTEGER const ldw, REAL *rwork, REAL &rcond, REAL &resid) {
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
    REAL anorm = Clansp("1", uplo, n, a, rwork);
    REAL ainvnm = Clansp("1", uplo, n, ainv, rwork);
    if (anorm <= zero || ainvnm <= zero) {
        rcond = zero;
        resid = one / eps;
        return;
    }
    rcond = (one / anorm) / ainvnm;
    //
    //     Case where both A and AINV are upper triangular:
    //     Each element of - A * AINV is computed by taking the dot product
    //     of a row of A with a column of AINV.
    //
    INTEGER i = 0;
    INTEGER icol = 0;
    INTEGER j = 0;
    INTEGER jcol = 0;
    COMPLEX t = 0.0;
    INTEGER kcol = 0;
    INTEGER k = 0;
    INTEGER nall = 0;
    if (Mlsame(uplo, "U")) {
        for (i = 1; i <= n; i = i + 1) {
            icol = ((i - 1) * i) / 2 + 1;
            //
            //           Code when J <= I
            //
            for (j = 1; j <= i; j = j + 1) {
                jcol = ((j - 1) * j) / 2 + 1;
                t = Cdotu(j, &a[icol - 1], 1, ainv[jcol - 1], 1);
                jcol += 2 * j - 1;
                kcol = icol - 1;
                for (k = j + 1; k <= i; k = k + 1) {
                    t += a[(kcol + k) - 1] * ainv[jcol - 1];
                    jcol += k;
                }
                kcol += 2 * i;
                for (k = i + 1; k <= n; k = k + 1) {
                    t += a[kcol - 1] * ainv[jcol - 1];
                    kcol += k;
                    jcol += k;
                }
                work[(i - 1) + (j - 1) * ldwork] = -t;
            }
            //
            //           Code when J > I
            //
            for (j = i + 1; j <= n; j = j + 1) {
                jcol = ((j - 1) * j) / 2 + 1;
                t = Cdotu(i, &a[icol - 1], 1, ainv[jcol - 1], 1);
                jcol = jcol - 1;
                kcol = icol + 2 * i - 1;
                for (k = i + 1; k <= j; k = k + 1) {
                    t += a[kcol - 1] * ainv[(jcol + k) - 1];
                    kcol += k;
                }
                jcol += 2 * j;
                for (k = j + 1; k <= n; k = k + 1) {
                    t += a[kcol - 1] * ainv[jcol - 1];
                    kcol += k;
                    jcol += k;
                }
                work[(i - 1) + (j - 1) * ldwork] = -t;
            }
        }
    } else {
        //
        //        Case where both A and AINV are lower triangular
        //
        nall = (n * (n + 1)) / 2;
        for (i = 1; i <= n; i = i + 1) {
            //
            //           Code when J <= I
            //
            icol = nall - ((n - i + 1) * (n - i + 2)) / 2 + 1;
            for (j = 1; j <= i; j = j + 1) {
                jcol = nall - ((n - j) * (n - j + 1)) / 2 - (n - i);
                t = Cdotu(n - i + 1, &a[icol - 1], 1, ainv[jcol - 1], 1);
                kcol = i;
                jcol = j;
                for (k = 1; k <= j - 1; k = k + 1) {
                    t += a[kcol - 1] * ainv[jcol - 1];
                    jcol += n - k;
                    kcol += n - k;
                }
                jcol = jcol - j;
                for (k = j; k <= i - 1; k = k + 1) {
                    t += a[kcol - 1] * ainv[(jcol + k) - 1];
                    kcol += n - k;
                }
                work[(i - 1) + (j - 1) * ldwork] = -t;
            }
            //
            //           Code when J > I
            //
            icol = nall - ((n - i) * (n - i + 1)) / 2;
            for (j = i + 1; j <= n; j = j + 1) {
                jcol = nall - ((n - j + 1) * (n - j + 2)) / 2 + 1;
                t = Cdotu(n - j + 1, &a[(icol - n + j) - 1], 1, ainv[jcol - 1], 1);
                kcol = i;
                jcol = j;
                for (k = 1; k <= i - 1; k = k + 1) {
                    t += a[kcol - 1] * ainv[jcol - 1];
                    jcol += n - k;
                    kcol += n - k;
                }
                kcol = kcol - i;
                for (k = i; k <= j - 1; k = k + 1) {
                    t += a[(kcol + k) - 1] * ainv[jcol - 1];
                    jcol += n - k;
                }
                work[(i - 1) + (j - 1) * ldwork] = -t;
            }
        }
    }
    //
    //     Add the identity matrix to WORK .
    //
    for (i = 1; i <= n; i = i + 1) {
        work[(i - 1) + (i - 1) * ldwork] += one;
    }
    //
    //     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
    //
    resid = Clange("1", n, n, work, ldw, rwork);
    //
    resid = ((resid * rcond) / eps) / n.real();
    //
    //     End of Cspt03
    //
}
