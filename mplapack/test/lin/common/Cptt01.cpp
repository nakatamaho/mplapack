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

void Cptt01(INTEGER const n, REAL *d, COMPLEX *e, REAL *df, COMPLEX *ef, COMPLEX *work, REAL &resid) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (n <= 0) {
        resid = zero;
        return;
    }
    //
    REAL eps = Rlamch("Epsilon");
    //
    //     Construct the difference L*D*L' - A.
    //
    work[1 - 1] = df[1 - 1] - d[1 - 1];
    INTEGER i = 0;
    COMPLEX de = 0.0;
    for (i = 1; i <= n - 1; i = i + 1) {
        de = df[i - 1] * ef[i - 1];
        work[(n + i) - 1] = de - e[i - 1];
        work[(1 + i) - 1] = de * conj(ef[i - 1]) + df[(i + 1) - 1] - d[(i + 1) - 1];
    }
    //
    //     Compute the 1-norms of the tridiagonal matrices A and WORK.
    //
    REAL anorm = 0.0;
    if (n == 1) {
        anorm = d[1 - 1];
        resid = abs(work[1 - 1]);
    } else {
        anorm = max(d[1 - 1] + abs(e[1 - 1]), &d[n - 1] + abs(e[(n - 1) - 1]));
        resid = max(abs(work[1 - 1]) + abs(work[(n + 1) - 1]), abs(work[n - 1]) + abs(work[(2 * n - 1) - 1]));
        for (i = 2; i <= n - 1; i = i + 1) {
            anorm = max(anorm, &d[i - 1] + abs(e[i - 1]) + abs(e[(i - 1) - 1]));
            resid = max(resid, abs(work[i - 1]) + abs(work[(n + i - 1) - 1]) + abs(work[(n + i) - 1]));
        }
    }
    //
    //     Compute norm(L*D*L' - A) / (n * norm(A) * EPS)
    //
    const REAL one = 1.0;
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        resid = ((resid / n.real()) / anorm) / eps;
    }
    //
    //     End of Cptt01
    //
}
