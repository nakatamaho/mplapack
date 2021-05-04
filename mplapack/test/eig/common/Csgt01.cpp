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

void Csgt01(INTEGER const itype, const char *uplo, INTEGER const n, INTEGER const m, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *z, INTEGER const ldz, REAL *d, COMPLEX *work, REAL *rwork, REAL *result) {
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
    REAL anorm = Clanhe("1", uplo, n, a, lda, rwork) * Clange("1", n, m, z, ldz, rwork);
    const REAL one = 1.0;
    if (anorm == zero) {
        anorm = one;
    }
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER i = 0;
    if (itype == 1) {
        //
        //        Norm of AZ - BZD
        //
        Chemm("Left", uplo, n, m, cone, a, lda, z, ldz, czero, work, n);
        for (i = 1; i <= m; i = i + 1) {
            CRscal(n, &d[i - 1], &z[(i - 1) * ldz], 1);
        }
        Chemm("Left", uplo, n, m, cone, b, ldb, z, ldz, -cone, work, n);
        //
        result[1 - 1] = (Clange("1", n, m, work, n, rwork) / anorm) / (n * ulp);
        //
    } else if (itype == 2) {
        //
        //        Norm of ABZ - ZD
        //
        Chemm("Left", uplo, n, m, cone, b, ldb, z, ldz, czero, work, n);
        for (i = 1; i <= m; i = i + 1) {
            CRscal(n, &d[i - 1], &z[(i - 1) * ldz], 1);
        }
        Chemm("Left", uplo, n, m, cone, a, lda, work, n, -cone, z, ldz);
        //
        result[1 - 1] = (Clange("1", n, m, z, ldz, rwork) / anorm) / (n * ulp);
        //
    } else if (itype == 3) {
        //
        //        Norm of BAZ - ZD
        //
        Chemm("Left", uplo, n, m, cone, a, lda, z, ldz, czero, work, n);
        for (i = 1; i <= m; i = i + 1) {
            CRscal(n, &d[i - 1], &z[(i - 1) * ldz], 1);
        }
        Chemm("Left", uplo, n, m, cone, b, ldb, work, n, -cone, z, ldz);
        //
        result[1 - 1] = (Clange("1", n, m, z, ldz, rwork) / anorm) / (n * ulp);
    }
    //
    //     End of CDGT01
    //
}
