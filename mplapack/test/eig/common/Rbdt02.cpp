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

void Rbdt02(INTEGER const m, INTEGER const n, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL *u, INTEGER const ldu, REAL *work, REAL &resid) {
    b([ldb * star]);
    c([ldc * star]);
    u([ldu * star]);
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
    // ======================================================================
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
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    resid = zero;
    if (m <= 0 || n <= 0) {
        return;
    }
    REAL realmn = (max(m, n)).real();
    REAL eps = Rlamch("Precision");
    //
    //     Compute norm( B - U * C )
    //
    INTEGER j = 0;
    const REAL one = 1.0;
    for (j = 1; j <= n; j = j + 1) {
        Rcopy(m, &b[(j - 1) * ldb], 1, work, 1);
        Rgemv("No transpose", m, m, -one, u, ldu, &c[(j - 1) * ldc], 1, one, work, 1);
        resid = max({resid, Rasum(m, work, 1)});
    }
    //
    //     Compute norm of B.
    //
    REAL bnorm = Rlange("1", m, n, b, ldb, work);
    //
    if (bnorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        if (bnorm >= resid) {
            resid = (resid / bnorm) / (realmn * eps);
        } else {
            if (bnorm < one) {
                resid = (min(resid, realmn * bnorm) / bnorm) / (realmn * eps);
            } else {
                resid = min(resid / bnorm, realmn) / (realmn * eps);
            }
        }
    }
    //
    //     End of Rbdt02
    //
}
