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

void Rptts2(INTEGER const n, INTEGER const nrhs, REAL *d, REAL *e, REAL *b, INTEGER const ldb) {
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
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        if (n == 1) {
            Rscal(nrhs, 1.0 / d[1 - 1], b, ldb);
        }
        return;
    }
    //
    //     Solve A * X = B using the factorization A = L*D*L**T,
    //     overwriting each right hand side vector with its solution.
    //
    INTEGER j = 0;
    INTEGER i = 0;
    for (j = 1; j <= nrhs; j = j + 1) {
        //
        //           Solve L * x = b.
        //
        for (i = 2; i <= n; i = i + 1) {
            b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - b[((i - 1) - 1) + (j - 1) * ldb] * e[(i - 1) - 1];
        }
        //
        //           Solve D * L**T * x = b.
        //
        b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
        for (i = n - 1; i >= 1; i = i - 1) {
            b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] / d[i - 1] - b[((i + 1) - 1) + (j - 1) * ldb] * e[i - 1];
        }
    }
    //
    //     End of Rptts2
    //
}
