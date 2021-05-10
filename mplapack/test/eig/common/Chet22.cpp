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

void Chet22(INTEGER const itype, const char *uplo, INTEGER const n, INTEGER const m, INTEGER const kband, COMPLEX *a, INTEGER const lda, REAL *d, REAL *e, COMPLEX *u, INTEGER const ldu, COMPLEX * /* v */, INTEGER const ldv, COMPLEX * /* tau */, COMPLEX *work, REAL *rwork, REAL *result) {
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
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    result[2 - 1] = zero;
    if (n <= 0 || m <= 0) {
        return;
    }
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Precision");
    //
    //     Do Test 1
    //
    //     Norm of A:
    //
    REAL anorm = max({Clanhe("1", uplo, n, a, lda, rwork), unfl});
    //
    //     Compute error matrix:
    //
    //     ITYPE=1: error = U**H A U - S
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Chemm("L", uplo, n, m, cone, a, lda, u, ldu, czero, work, n);
    INTEGER nn = n * n;
    INTEGER nnp1 = nn + 1;
    Cgemm("C", "N", m, m, n, cone, u, ldu, work, n, czero, &work[nnp1 - 1], n);
    INTEGER j = 0;
    INTEGER jj = 0;
    for (j = 1; j <= m; j = j + 1) {
        jj = nn + (j - 1) * n + j;
        work[jj - 1] = work[jj - 1] - d[j - 1];
    }
    INTEGER jj1 = 0;
    INTEGER jj2 = 0;
    if (kband == 1 && n > 1) {
        for (j = 2; j <= m; j = j + 1) {
            jj1 = nn + (j - 1) * n + j - 1;
            jj2 = nn + (j - 2) * n + j;
            work[jj1 - 1] = work[jj1 - 1] - e[(j - 1) - 1];
            work[jj2 - 1] = work[jj2 - 1] - e[(j - 1) - 1];
        }
    }
    REAL wnorm = Clanhe("1", uplo, m, &work[nnp1 - 1], n, rwork);
    //
    const REAL one = 1.0;
    if (anorm > wnorm) {
        result[1 - 1] = (wnorm / anorm) / (m * ulp);
    } else {
        if (anorm < one) {
            result[1 - 1] = (min(wnorm, m * anorm) / anorm) / (m * ulp);
        } else {
            result[1 - 1] = min(wnorm / anorm, castREAL(m)) / (m * ulp);
        }
    }
    //
    //     Do Test 2
    //
    //     Compute  U**H U - I
    //
    if (itype == 1) {
        Cunt01("Columns", n, m, u, ldu, work, 2 * n * n, rwork, result[2 - 1]);
    }
    //
    //     End of Chet22
    //
}
