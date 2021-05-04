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

void Rget51(INTEGER const itype, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *u, INTEGER const ldu, REAL *v, INTEGER const ldv, REAL *work, REAL &result) {
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
    const REAL one = 1.0;
    INTEGER jcol = 0;
    INTEGER jrow = 0;
    REAL wnorm = 0.0;
    INTEGER jdiag = 0;
    if (itype <= 2) {
        //
        //        Tests scaled by the norm(A)
        //
        anorm = max({Rlange("1", n, n, a, lda, work), unfl});
        //
        if (itype == 1) {
            //
            //           ITYPE=1: Compute W = A - UBV'
            //
            Rlacpy(" ", n, n, a, lda, work, n);
            Rgemm("N", "N", n, n, n, one, u, ldu, b, ldb, zero, &work[(pow2(n) + 1) - 1], n);
            //
            Rgemm("N", "C", n, n, n, -one, &work[(pow2(n) + 1) - 1], n, v, ldv, one, work, n);
            //
        } else {
            //
            //           ITYPE=2: Compute W = A - B
            //
            Rlacpy(" ", n, n, b, ldb, work, n);
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
        wnorm = Rlange("1", n, n, work, n, &work[(pow2(n) + 1) - 1]);
        //
        if (anorm > wnorm) {
            result = (wnorm / anorm) / (n * ulp);
        } else {
            if (anorm < one) {
                result = (min(wnorm, n * anorm) / anorm) / (n * ulp);
            } else {
                result = min(wnorm / anorm, n.real()) / (n * ulp);
            }
        }
        //
    } else {
        //
        //        Tests not scaled by norm(A)
        //
        //        ITYPE=3: Compute  UU' - I
        //
        Rgemm("N", "C", n, n, n, one, u, ldu, u, ldu, zero, work, n);
        //
        for (jdiag = 1; jdiag <= n; jdiag = jdiag + 1) {
            work[((n + 1) * (jdiag - 1) + 1) - 1] = work[((n + 1) * (jdiag - 1) + 1) - 1] - one;
        }
        //
        result = min({Rlange("1", n, n, work, n, &work[(pow2(n) + 1) - 1]), n.real()}) / (n * ulp);
    }
    //
    //     End of Rget51
    //
}
