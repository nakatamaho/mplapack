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

void Rbdt05(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *s, INTEGER const ns, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, REAL *work, REAL &resid) {
    a([lda * star]);
    u([ldu * star]);
    vt([ldvt * star]);
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
    //     Quick return if possible.
    //
    const REAL zero = 0.0;
    resid = zero;
    if (min(m, n) <= 0 || ns <= 0) {
        return;
    }
    //
    REAL eps = Rlamch("Precision");
    REAL anorm = Rlange("M", m, n, a, lda, work);
    //
    //     Compute U' * A * V.
    //
    const REAL one = 1.0;
    Rgemm("N", "T", m, ns, n, one, a, lda, vt, ldvt, zero, &work[(1 + ns * ns) - 1], m);
    Rgemm("T", "N", ns, ns, m, -one, u, ldu, &work[(1 + ns * ns) - 1], m, zero, work, ns);
    //
    //     norm(S - U' * B * V)
    //
    INTEGER j = 0;
    INTEGER i = 0;
    for (i = 1; i <= ns; i = i + 1) {
        work[(j + i) - 1] += s[i - 1];
        resid = max({resid, Rasum(ns, &work[(j + 1) - 1], 1)});
        j += ns;
    }
    //
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        if (anorm >= resid) {
            resid = (resid / anorm) / (n.real() * eps);
        } else {
            if (anorm < one) {
                resid = (min(resid, n.real() * anorm) / anorm) / (n.real() * eps);
            } else {
                resid = min(resid / anorm, n.real()) / (n.real() * eps);
            }
        }
    }
    //
    //     End of Rbdt05
    //
}
