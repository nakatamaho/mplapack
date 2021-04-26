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

void Rlarge(INTEGER const n, REAL *a, INTEGER const lda, INTEGER *iseed, REAL *work, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (lda < max((INTEGER)1, n)) {
        info = -3;
    }
    if (info < 0) {
        Mxerbla("Rlarge", -info);
        return;
    }
    //
    //     pre- and post-multiply A by random orthogonal matrix
    //
    INTEGER i = 0;
    REAL wn = 0.0;
    REAL wa = 0.0;
    const REAL zero = 0.0;
    REAL tau = 0.0;
    REAL wb = 0.0;
    const REAL one = 1.0;
    for (i = n; i >= 1; i = i - 1) {
        //
        //        generate random reflection
        //
        Rlarnv(3, iseed, n - i + 1, work);
        wn = Rnrm2(n - i + 1, work, 1);
        wa = sign(wn, work[1 - 1]);
        if (wn == zero) {
            tau = zero;
        } else {
            wb = work[1 - 1] + wa;
            Rscal(n - i, one / wb, &work[2 - 1], 1);
            work[1 - 1] = one;
            tau = wb / wa;
        }
        //
        //        multiply A(i:n,1:n) by random reflection from the left
        //
        Rgemv("Transpose", n - i + 1, n, one, &a[(i - 1)], lda, work, 1, zero, &work[(n + 1) - 1], 1);
        Rger(n - i + 1, n, -tau, work, 1, &work[(n + 1) - 1], 1, &a[(i - 1)], lda);
        //
        //        multiply A(1:n,i:n) by random reflection from the right
        //
        Rgemv("No transpose", n, n - i + 1, one, &a[(i - 1) * lda], lda, work, 1, zero, &work[(n + 1) - 1], 1);
        Rger(n, n - i + 1, -tau, &work[(n + 1) - 1], 1, work, 1, &a[(i - 1) * lda], lda);
    }
    //
    //     End of Rlarge
    //
}
