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

void Cbdt02(INTEGER const m, INTEGER const n, COMPLEX *b, INTEGER const ldb, COMPLEX *c, INTEGER const ldc, COMPLEX *u, INTEGER const ldu, COMPLEX *work, REAL *rwork, REAL &resid) {
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
    REAL realmn = castREAL(max(m, n));
    REAL eps = Rlamch("Precision");
    //
    //     Compute norm( B - U * C )
    //
    INTEGER j = 0;
    const REAL one = 1.0;
    for (j = 1; j <= n; j = j + 1) {
        Ccopy(m, &b[(j - 1) * ldb], 1, work, 1);
        Cgemv("No transpose", m, m, -COMPLEX(one), u, ldu, &c[(j - 1) * ldc], 1, COMPLEX(one), work, 1);
        resid = max({resid, RCasum(m, work, 1)});
    }
    //
    //     Compute norm of B.
    //
    REAL bnorm = Clange("1", m, n, b, ldb, rwork);
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
    //     End of Cbdt02
    //
}
