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
#include <mplapack.h>

void Clarcm(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *c, INTEGER const ldc, REAL *rwork) {
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible.
    //
    if ((m == 0) || (n == 0)) {
        return;
    }
    //
    INTEGER j = 0;
    INTEGER i = 0;
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            rwork[((j - 1) * m + i) - 1] = b[(i - 1) + (j - 1) * ldb].real();
        }
    }
    //
    INTEGER l = m * n + 1;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    Rgemm("N", "N", m, n, m, one, a, lda, rwork, m, zero, &rwork[l - 1], m);
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            c[(i - 1) + (j - 1) * ldc] = rwork[(l + (j - 1) * m + i - 1) - 1];
        }
    }
    //
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            rwork[((j - 1) * m + i) - 1] = b[(i - 1) + (j - 1) * ldb].imag();
        }
    }
    Rgemm("N", "N", m, n, m, one, a, lda, rwork, m, zero, &rwork[l - 1], m);
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            c[(i - 1) + (j - 1) * ldc] = COMPLEX(c[(i - 1) + (j - 1) * ldc].real(), &rwork[(l + (j - 1) * m + i - 1) - 1]);
        }
    }
    //
    //     End of Clarcm
    //
}
