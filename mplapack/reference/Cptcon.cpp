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

void Cptcon(INTEGER const &n, REAL *d, COMPLEX *e, REAL const &anorm, REAL &rcond, REAL *rwork, INTEGER &info) {
    //
    //  -- LAPACK computational routine --
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
    //     .. Executable Statements ..
    //
    //     Test the input arguments.
    //
    info = 0;
    const REAL zero = 0.0;
    if (n < 0) {
        info = -1;
    } else if (anorm < zero) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Cptcon", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    rcond = zero;
    const REAL one = 1.0;
    if (n == 0) {
        rcond = one;
        return;
    } else if (anorm == zero) {
        return;
    }
    //
    //     Check that D(1:N) is positive.
    //
    INTEGER i = 0;
    for (i = 1; i <= n; i = i + 1) {
        if (d[i - 1] <= zero) {
            return;
        }
    }
    //
    //     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
    //
    //        m(i,j) =  abs(A(i,j)), i = j,
    //        m(i,j) = -abs(A(i,j)), i .ne. j,
    //
    //     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.
    //
    //     Solve M(L) * x = e.
    //
    rwork[1 - 1] = one;
    for (i = 2; i <= n; i = i + 1) {
        rwork[i - 1] = one + rwork[(i - 1) - 1] * abs(e[(i - 1) - 1]);
    }
    //
    //     Solve D * M(L)**H * x = b.
    //
    rwork[n - 1] = rwork[n - 1] / d[n - 1];
    for (i = n - 1; i >= 1; i = i - 1) {
        rwork[i - 1] = rwork[i - 1] / d[i - 1] + rwork[(i + 1) - 1] * abs(e[i - 1]);
    }
    //
    //     Compute AINVNM = max(x(i)), 1<=i<=n.
    //
    INTEGER ix = iRamax[(n - 1) + (rwork - 1) * ldiRamax];
    REAL ainvnm = abs(rwork[ix - 1]);
    //
    //     Compute the reciprocal condition number.
    //
    if (ainvnm != zero) {
        rcond = (one / ainvnm) / anorm;
    }
    //
    //     End of Cptcon
    //
}
