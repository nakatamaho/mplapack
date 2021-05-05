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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Rsgt01(INTEGER const itype, const char *uplo, INTEGER const n, INTEGER const m, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *z, INTEGER const ldz, REAL *d, REAL *work, REAL *result) {
    a([lda * star]);
    b([ldb * star]);
    z([ldz * star]);
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
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    if (n <= 0) {
        return;
    }
    //
    REAL ulp = Rlamch("Epsilon");
    //
    //     Compute product of 1-norms of A and Z.
    //
    REAL anorm = Rlansy("1", uplo, n, a, lda, work) * Rlange("1", n, m, z, ldz, work);
    const REAL one = 1.0;
    if (anorm == zero) {
        anorm = one;
    }
    //
    INTEGER i = 0;
    if (itype == 1) {
        //
        //        Norm of AZ - BZD
        //
        Rsymm("Left", uplo, n, m, one, a, lda, z, ldz, zero, work, n);
        for (i = 1; i <= m; i = i + 1) {
            Rscal(n, &d[i - 1], &z[(i - 1) * ldz], 1);
        }
        Rsymm("Left", uplo, n, m, one, b, ldb, z, ldz, -one, work, n);
        //
        result[1 - 1] = (Rlange("1", n, m, work, n, work) / anorm) / (n * ulp);
        //
    } else if (itype == 2) {
        //
        //        Norm of ABZ - ZD
        //
        Rsymm("Left", uplo, n, m, one, b, ldb, z, ldz, zero, work, n);
        for (i = 1; i <= m; i = i + 1) {
            Rscal(n, &d[i - 1], &z[(i - 1) * ldz], 1);
        }
        Rsymm("Left", uplo, n, m, one, a, lda, work, n, -one, z, ldz);
        //
        result[1 - 1] = (Rlange("1", n, m, z, ldz, work) / anorm) / (n * ulp);
        //
    } else if (itype == 3) {
        //
        //        Norm of BAZ - ZD
        //
        Rsymm("Left", uplo, n, m, one, a, lda, z, ldz, zero, work, n);
        for (i = 1; i <= m; i = i + 1) {
            Rscal(n, &d[i - 1], &z[(i - 1) * ldz], 1);
        }
        Rsymm("Left", uplo, n, m, one, b, ldb, work, n, -one, z, ldz);
        //
        result[1 - 1] = (Rlange("1", n, m, z, ldz, work) / anorm) / (n * ulp);
    }
    //
    //     End of DDGT01
    //
}
