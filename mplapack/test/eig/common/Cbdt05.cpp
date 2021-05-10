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

void Cbdt05(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *s, INTEGER const ns, COMPLEX *u, INTEGER const ldu, COMPLEX *vt, INTEGER const ldvt, COMPLEX *work, REAL &resid) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible.
    //
    const REAL zero = 0.0;
    resid = zero;
    if (min(m, n) <= 0 || ns <= 0) {
        return;
    }
    //
    REAL eps = Rlamch("Precision");
    REAL dum[1];
    REAL anorm = Clange("M", m, n, a, lda, dum);
    //
    //     Compute U' * A * V.
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Cgemm("N", "C", m, ns, n, cone, a, lda, vt, ldvt, czero, &work[(1 + ns * ns) - 1], m);
    Cgemm("C", "N", ns, ns, m, -cone, u, ldu, &work[(1 + ns * ns) - 1], m, czero, work, ns);
    //
    //     norm(S - U' * B * V)
    //
    INTEGER j = 0;
    INTEGER i = 0;
    for (i = 1; i <= ns; i = i + 1) {
        work[(j + i) - 1] += COMPLEX(s[i - 1], zero);
        resid = max({resid, RCasum(ns, &work[(j + 1) - 1], 1)});
        j += ns;
    }
    //
    const REAL one = 1.0;
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        if (anorm >= resid) {
            resid = (resid / anorm) / (castREAL(n) * eps);
        } else {
            if (anorm < one) {
                resid = (min(resid, castREAL(n) * anorm) / anorm) / (castREAL(n) * eps);
            } else {
                resid = min(resid / anorm, castREAL(n)) / (castREAL(n) * eps);
            }
        }
    }
    //
    //     End of Cbdt05
    //
}
