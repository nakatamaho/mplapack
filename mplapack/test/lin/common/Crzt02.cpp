/*
 * Copyright (c) 2008-2021
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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

REAL Crzt02(INTEGER const m, INTEGER const n, COMPLEX *af, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork) {
    REAL return_value = 0.0;
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
    const REAL zero = 0.0;
    return_value = zero;
    //
    if (lwork < n * n + n) {
        Mxerbla("Crzt02", 7);
        return return_value;
    }
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0) {
        return return_value;
    }
    //
    //     Q := I
    //
    const REAL one = 1.0;
    Claset("Full", n, n, COMPLEX(zero), COMPLEX(one), work, n);
    //
    //     Q := P(1) * ... * P(m) * Q
    //
    INTEGER info = 0;
    Cunmrz("Left", "No transpose", n, n, m, n - m, af, lda, tau, work, n, &work[(n * n + 1) - 1], lwork - n * n, info);
    //
    //     Q := P(m)' * ... * P(1)' * Q
    //
    Cunmrz("Left", "Conjugate transpose", n, n, m, n - m, af, lda, tau, work, n, &work[(n * n + 1) - 1], lwork - n * n, info);
    //
    //     Q := Q - I
    //
    INTEGER i = 0;
    for (i = 1; i <= n; i = i + 1) {
        work[((i - 1) * n + i) - 1] = work[((i - 1) * n + i) - 1] - one;
    }
    //
    REAL rwork[1];
    return_value = Clange("One-norm", n, n, work, n, rwork) / (Rlamch("Epsilon") * (max(m, n)).real());
    return return_value;
    //
    //     End of Crzt02
    //
}
