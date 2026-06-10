/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZTPT01.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Ctpt01(fem::str_cref uplo, fem::str_cref diag, INTEGER const n, COMPLEX *ap, COMPLEX *ainvp, REAL &rcond, REAL *rwork, REAL &resid) {
    //
    // Quick exit if N = 0.
    //
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (n <= 0) {
        rcond = one;
        resid = zero;
        return;
    }
    //
    // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
    //
    REAL eps = Rlamch("Epsilon");
    REAL anorm = Clantp("1", uplo.elems(), diag.elems(), n, ap, rwork);
    REAL ainvnm = Clantp("1", uplo.elems(), diag.elems(), n, ainvp, rwork);
    if (anorm <= zero || ainvnm <= zero) {
        rcond = zero;
        resid = one / eps;
        return;
    }
    rcond = (one / anorm) / ainvnm;
    //
    // Compute A * AINV, overwriting AINV.
    //
    bool unitd = Mlsame(diag.elems(), "U");
    INTEGER jc = 0;
    INTEGER j = 0;
    if (Mlsame(uplo.elems(), "U")) {
        jc = 1;
        for (j = 1; j <= n; j = j + 1) {
            if (unitd) {
                ainvp[(jc + j - 1) - 1] = one;
            }
            //
            // Form the j-th column of A*AINV.
            //
            Ctpmv("Upper", "No transpose", diag.elems(), j, ap, &ainvp[jc - 1], 1);
            //
            // Subtract 1 from the diagonal to form A*AINV - I.
            //
            ainvp[(jc + j - 1) - 1] = ainvp[(jc + j - 1) - 1] - one;
            jc += j;
        }
    } else {
        jc = 1;
        for (j = 1; j <= n; j = j + 1) {
            if (unitd) {
                ainvp[jc - 1] = one;
            }
            //
            // Form the j-th column of A*AINV.
            //
            Ctpmv("Lower", "No transpose", diag.elems(), n - j + 1, &ap[jc - 1], &ainvp[jc - 1], 1);
            //
            // Subtract 1 from the diagonal to form A*AINV - I.
            //
            ainvp[jc - 1] = ainvp[jc - 1] - one;
            jc += n - j + 1;
        }
    }
    //
    // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
    //
    resid = Clantp("1", uplo.elems(), "Non-unit", n, ainvp, rwork);
    //
    resid = ((resid * rcond) / castREAL(n)) / eps;
    //
    // End of Ctpt01
    //
}
