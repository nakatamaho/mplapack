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

void Rgtts2(INTEGER const &itrans, INTEGER const &n, INTEGER const &nrhs, REAL *dl, REAL *d, REAL *du, REAL *du2, INTEGER *ipiv, REAL *b, INTEGER const &ldb) {
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER ip = 0;
    REAL temp = 0.0;
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
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    if (itrans == 0) {
        //
        //        Solve A*X = B using the LU factorization of A,
        //        overwriting each right hand side vector with its solution.
        //
        if (nrhs <= 1) {
            j = 1;
        statement_10:
            //
            //           Solve L*x = b.
            //
            for (i = 1; i <= n - 1; i = i + 1) {
                ip = ipiv[i - 1];
                temp = b[((i + 1 - ip + i) - 1) + (j - 1) * ldb] - dl[i - 1] * b[(ip - 1) + (j - 1) * ldb];
                b[(i - 1) + (j - 1) * ldb] = b[(ip - 1) + (j - 1) * ldb];
                b[((i + 1) - 1) + (j - 1) * ldb] = temp;
            }
            //
            //           Solve U*x = b.
            //
            b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
            if (n > 1) {
                b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
            }
            for (i = n - 2; i >= 1; i = i - 1) {
                b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb] - du2[i - 1] * b[((i + 2) - 1) + (j - 1) * ldb]) / d[i - 1];
            }
            if (j < nrhs) {
                j++;
                goto statement_10;
            }
        } else {
            for (j = 1; j <= nrhs; j = j + 1) {
                //
                //              Solve L*x = b.
                //
                for (i = 1; i <= n - 1; i = i + 1) {
                    if (ipiv[i - 1] == i) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - dl[i - 1] * b[(i - 1) + (j - 1) * ldb];
                    } else {
                        temp = b[(i - 1) + (j - 1) * ldb];
                        b[(i - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb];
                        b[((i + 1) - 1) + (j - 1) * ldb] = temp - dl[i - 1] * b[(i - 1) + (j - 1) * ldb];
                    }
                }
                //
                //              Solve U*x = b.
                //
                b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
                if (n > 1) {
                    b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
                }
                for (i = n - 2; i >= 1; i = i - 1) {
                    b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb] - du2[i - 1] * b[((i + 2) - 1) + (j - 1) * ldb]) / d[i - 1];
                }
            }
        }
    } else {
        //
        //        Solve A**T * X = B.
        //
        if (nrhs <= 1) {
            //
            //           Solve U**T*x = b.
            //
            j = 1;
        statement_70:
            b[(j - 1) * ldb] = b[(j - 1) * ldb] / d[1 - 1];
            if (n > 1) {
                b[(2 - 1) + (j - 1) * ldb] = (b[(2 - 1) + (j - 1) * ldb] - du[1 - 1] * b[(j - 1) * ldb]) / d[2 - 1];
            }
            for (i = 3; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[(i - 1) - 1] * b[((i - 1) - 1) + (j - 1) * ldb] - du2[(i - 2) - 1] * b[((i - 2) - 1) + (j - 1) * ldb]) / d[i - 1];
            }
            //
            //           Solve L**T*x = b.
            //
            for (i = n - 1; i >= 1; i = i - 1) {
                ip = ipiv[i - 1];
                temp = b[(i - 1) + (j - 1) * ldb] - dl[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb];
                b[(i - 1) + (j - 1) * ldb] = b[(ip - 1) + (j - 1) * ldb];
                b[(ip - 1) + (j - 1) * ldb] = temp;
            }
            if (j < nrhs) {
                j++;
                goto statement_70;
            }
            //
        } else {
            for (j = 1; j <= nrhs; j = j + 1) {
                //
                //              Solve U**T*x = b.
                //
                b[(j - 1) * ldb] = b[(j - 1) * ldb] / d[1 - 1];
                if (n > 1) {
                    b[(2 - 1) + (j - 1) * ldb] = (b[(2 - 1) + (j - 1) * ldb] - du[1 - 1] * b[(j - 1) * ldb]) / d[2 - 1];
                }
                for (i = 3; i <= n; i = i + 1) {
                    b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[(i - 1) - 1] * b[((i - 1) - 1) + (j - 1) * ldb] - du2[(i - 2) - 1] * b[((i - 2) - 1) + (j - 1) * ldb]) / d[i - 1];
                }
                for (i = n - 1; i >= 1; i = i - 1) {
                    if (ipiv[i - 1] == i) {
                        b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - dl[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb];
                    } else {
                        temp = b[((i + 1) - 1) + (j - 1) * ldb];
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - dl[i - 1] * temp;
                        b[(i - 1) + (j - 1) * ldb] = temp;
                    }
                }
            }
        }
    }
    //
    //     End of Rgtts2
    //
}
