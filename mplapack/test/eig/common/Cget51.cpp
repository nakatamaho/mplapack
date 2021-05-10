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

void Cget51(INTEGER const itype, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, COMPLEX *work, REAL *rwork, REAL &result) {
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
    result = zero;
    if (n <= 0) {
        return;
    }
    //
    //     Constants
    //
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    //
    //     Some Error Checks
    //
    const REAL ten = 10.0;
    if (itype < 1 || itype > 3) {
        result = ten / ulp;
        return;
    }
    //
    REAL anorm = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER jcol = 0;
    INTEGER jrow = 0;
    REAL wnorm = 0.0;
    const REAL one = 1.0;
    INTEGER jdiag = 0;
    if (itype <= 2) {
        //
        //        Tests scaled by the norm(A)
        //
        anorm = max({Clange("1", n, n, a, lda, rwork), unfl});
        //
        if (itype == 1) {
            //
            //           ITYPE=1: Compute W = A - U B V**H
            //
            Clacpy(" ", n, n, a, lda, work, n);
            Cgemm("N", "N", n, n, n, cone, u, ldu, b, ldb, czero, &work[(pow2(n) + 1) - 1], n);
            //
            Cgemm("N", "C", n, n, n, -cone, &work[(pow2(n) + 1) - 1], n, v, ldv, cone, work, n);
            //
        } else {
            //
            //           ITYPE=2: Compute W = A - B
            //
            Clacpy(" ", n, n, b, ldb, work, n);
            //
            for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                for (jrow = 1; jrow <= n; jrow = jrow + 1) {
                    work[(jrow + n * (jcol - 1)) - 1] = work[(jrow + n * (jcol - 1)) - 1] - a[(jrow - 1) + (jcol - 1) * lda];
                }
            }
        }
        //
        //        Compute norm(W)/ ( ulp*norm(A) )
        //
        wnorm = Clange("1", n, n, work, n, rwork);
        //
        if (anorm > wnorm) {
            result = (wnorm / anorm) / (n * ulp);
        } else {
            if (anorm < one) {
                result = (min(wnorm, n * anorm) / anorm) / (n * ulp);
            } else {
                result = min(wnorm / anorm, castREAL(n)) / (n * ulp);
            }
        }
        //
    } else {
        //
        //        Tests not scaled by norm(A)
        //
        //        ITYPE=3: Compute  U U**H - I
        //
        Cgemm("N", "C", n, n, n, cone, u, ldu, u, ldu, czero, work, n);
        //
        for (jdiag = 1; jdiag <= n; jdiag = jdiag + 1) {
            work[((n + 1) * (jdiag - 1) + 1) - 1] = work[((n + 1) * (jdiag - 1) + 1) - 1] - cone;
        }
        //
        result = min({Clange("1", n, n, work, n, rwork), castREAL(n)}) / (n * ulp);
    }
    //
    //     End of Cget51
    //
}
