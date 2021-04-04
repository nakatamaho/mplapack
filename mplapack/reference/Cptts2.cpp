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

void Cptts2(INTEGER const &iuplo, INTEGER const &n, INTEGER const &nrhs, REAL *d, COMPLEX *e, COMPLEX *b, INTEGER const &ldb) {
    INTEGER j = 0;
    INTEGER i = 0;
    //
    //  -- LAPACK computational routine --
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
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        if (n == 1) {
            CRscal(nrhs, 1.0 / d[1 - 1], b, ldb);
        }
        return;
    }
    //
    if (iuplo == 1) {
        //
        //        Solve A * X = B using the factorization A = U**H *D*U,
        //        overwriting each right hand side vector with its solution.
        //
        if (nrhs <= 2) {
            j = 1;
        statement_10:
            //
            //           Solve U**H * x = b.
            //
            for (i = 2; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i - 1) - 1) + (j - 1) * ldb] * conj(e[(i - 1) - 1]);
            }
            //
            //           Solve D * U * x = b.
            //
            for (i = 1; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] / d[i - 1];
            }
            for (i = n - 1; i >= 1; i = i - 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i + 1) - 1) + (j - 1) * ldb] * e[i - 1];
            }
            if (j < nrhs) {
                j++;
                goto statement_10;
            }
        } else {
            for (j = 1; j <= nrhs; j = j + 1) {
                //
                //              Solve U**H * x = b.
                //
                for (i = 2; i <= n; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i - 1) - 1) + (j - 1) * ldb] * conj(e[(i - 1) - 1]);
                }
                //
                //              Solve D * U * x = b.
                //
                b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
                for (i = n - 1; i >= 1; i = i - 1) {
                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] / d[i - 1] - b[((i + 1) - 1) + (j - 1) * ldb] * e[i - 1];
                }
            }
        }
    } else {
        //
        //        Solve A * X = B using the factorization A = L*D*L**H,
        //        overwriting each right hand side vector with its solution.
        //
        if (nrhs <= 2) {
            j = 1;
        statement_80:
            //
            //           Solve L * x = b.
            //
            for (i = 2; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i - 1) - 1) + (j - 1) * ldb] * e[(i - 1) - 1];
            }
            //
            //           Solve D * L**H * x = b.
            //
            for (i = 1; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] / d[i - 1];
            }
            for (i = n - 1; i >= 1; i = i - 1) {
                b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i + 1) - 1) + (j - 1) * ldb] * conj(e[i - 1]);
            }
            if (j < nrhs) {
                j++;
                goto statement_80;
            }
        } else {
            for (j = 1; j <= nrhs; j = j + 1) {
                //
                //              Solve L * x = b.
                //
                for (i = 2; i <= n; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i - 1) - 1) + (j - 1) * ldb] * e[(i - 1) - 1];
                }
                //
                //              Solve D * L**H * x = b.
                //
                b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
                for (i = n - 1; i >= 1; i = i - 1) {
                    b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] / d[i - 1] - b[((i + 1) - 1) + (j - 1) * ldb] * conj(e[i - 1]);
                }
            }
        }
    }
    //
    //     End of Cptts2
    //
}
