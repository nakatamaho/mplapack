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

void Rgesc2(INTEGER const &n, REAL *a, INTEGER const &lda, REAL *rhs, INTEGER *ipiv, INTEGER *jpiv, REAL &scale) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //      Set constant to control overflow
    //
    REAL eps = dlamch("P");
    REAL smlnum = dlamch("S") / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Apply permutations IPIV to RHS
    //
    Rlaswp(1, rhs, lda, 1, n - 1, ipiv, 1);
    //
    //     Solve for L part
    //
    INTEGER i = 0;
    INTEGER j = 0;
    for (i = 1; i <= n - 1; i = i + 1) {
        for (j = i + 1; j <= n; j = j + 1) {
            rhs[j - 1] = rhs[j - 1] - a[(j - 1) + (i - 1) * lda] * rhs[i - 1];
        }
    }
    //
    //     Solve for U part
    //
    scale = one;
    //
    //     Check for scaling
    //
    i = iRamax[(n - 1) + (rhs - 1) * ldiRamax];
    const REAL two = 2.0e+0;
    REAL temp = 0.0;
    if (two * smlnum * abs(rhs[i - 1]) > abs(a[(n - 1) + (n - 1) * lda])) {
        temp = (one / two) / abs(rhs[i - 1]);
        Rscal(n, temp, rhs[1 - 1], 1);
        scale = scale * temp;
    }
    //
    for (i = n; i >= 1; i = i - 1) {
        temp = one / a[(i - 1) + (i - 1) * lda];
        rhs[i - 1] = rhs[i - 1] * temp;
        for (j = i + 1; j <= n; j = j + 1) {
            rhs[i - 1] = rhs[i - 1] - rhs[j - 1] * (a[(i - 1) + (j - 1) * lda] * temp);
        }
    }
    //
    //     Apply permutations JPIV to the solution (RHS)
    //
    Rlaswp(1, rhs, lda, 1, n - 1, jpiv, -1);
    //
    //     End of Rgesc2
    //
}
