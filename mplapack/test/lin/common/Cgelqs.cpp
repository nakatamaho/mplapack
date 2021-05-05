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

void Cgelqs(INTEGER const m, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *b, INTEGER const ldb, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
    a([lda * star]);
    b([ldb * star]);
    work([lwork]);
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || m > n) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    } else if (lwork < 1 || lwork < nrhs && m > 0 && n > 0) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Cgelqs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0 || m == 0) {
        return;
    }
    //
    //     Solve L*X = B(1:m,:)
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    Ctrsm("Left", "Lower", "No transpose", "Non-unit", m, nrhs, cone, a, lda, b, ldb);
    //
    //     Set B(m+1:n,:) to zero
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (m < n) {
        Claset("Full", n - m, nrhs, czero, czero, &b[((m + 1) - 1)], ldb);
    }
    //
    //     B := Q' * B
    //
    Cunmlq("Left", "Conjugate transpose", n, nrhs, m, a, lda, tau, b, ldb, work, lwork, info);
    //
    //     End of Cgelqs
    //
}
