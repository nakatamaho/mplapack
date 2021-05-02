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
#include <mplapack_lin.h>

REAL Cqrt11(INTEGER const m, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    return_value = zero;
    //
    //     Test for sufficient workspace
    //
    if (lwork < m * m + m) {
        Mxerbla("Cqrt11", 7);
        return return_value;
    }
    //
    //     Quick return if possible
    //
    if (m <= 0) {
        return return_value;
    }
    //
    const REAL one = 1.0;
    Claset("Full", m, m, COMPLEX(zero), COMPLEX(one), work, m);
    //
    //     Form Q
    //
    INTEGER info = 0;
    Cunm2r("Left", "No transpose", m, m, k, a, lda, tau, work, m, &work[(m * m + 1) - 1], info);
    //
    //     Form Q'*Q
    //
    Cunm2r("Left", "Conjugate transpose", m, m, k, a, lda, tau, work, m, &work[(m * m + 1) - 1], info);
    //
    INTEGER j = 0;
    for (j = 1; j <= m; j = j + 1) {
        work[((j - 1) * m + j) - 1] = work[((j - 1) * m + j) - 1] - one;
    }
    //
    arr_1d<1, REAL> rdummy;
    return_value = Clange("One-norm", m, m, work, m, rdummy) / (m.real() * Rlamch("Epsilon"));
    //
    return return_value;
    //
    //     End of Cqrt11
    //
}
