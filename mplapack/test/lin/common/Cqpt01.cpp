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

REAL Cqpt01(INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, COMPLEX *af, INTEGER const lda, COMPLEX *tau, INTEGER *jpvt, COMPLEX *work, INTEGER const lwork) {
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
    //     Test if there is enough workspace
    //
    if (lwork < m * n + n) {
        Mxerbla("Cqpt01", 10);
        return return_value;
    }
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0) {
        return return_value;
    }
    //
    arr_1d<1, REAL> rwork;
    REAL norma = Clange("One-norm", m, n, a, lda, rwork);
    //
    INTEGER j = 0;
    INTEGER i = 0;
    for (j = 1; j <= k; j = j + 1) {
        for (i = 1; i <= min(j, m); i = i + 1) {
            work[((j - 1) * m + i) - 1] = af[(i - 1) + (j - 1) * ldaf];
        }
        for (i = j + 1; i <= m; i = i + 1) {
            work[((j - 1) * m + i) - 1] = zero;
        }
    }
    for (j = k + 1; j <= n; j = j + 1) {
        Ccopy(m, af[(j - 1) * ldaf], 1, &work[((j - 1) * m + 1) - 1], 1);
    }
    //
    INTEGER info = 0;
    Cunmqr("Left", "No transpose", m, n, k, af, lda, tau, work, m, &work[(m * n + 1) - 1], lwork - m * n, info);
    //
    const REAL one = 1.0;
    for (j = 1; j <= n; j = j + 1) {
        //
        //        Compare i-th column of QR and jpvt(i)-th column of A
        //
        Caxpy(m, COMPLEX(-one), &a[(jpvt[j - 1] - 1) * lda], 1, &work[((j - 1) * m + 1) - 1], 1);
    }
    //
    return_value = Clange("One-norm", m, n, work, m, rwork) / ((max(m, n)).real() * Rlamch("Epsilon"));
    if (norma != zero) {
        return_value = return_value / norma;
    }
    //
    return return_value;
    //
    //     End of Cqpt01
    //
}
