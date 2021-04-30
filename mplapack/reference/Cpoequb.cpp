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

void Cpoequb(INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *s, REAL &scond, REAL &amax, INTEGER &info) {
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
    //     Test the input parameters.
    //
    //     Positive definite only performs 1 pass of equilibration.
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (lda < max((INTEGER)1, n)) {
        info = -3;
    }
    if (info != 0) {
        Mxerbla("Cpoequb", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (n == 0) {
        scond = one;
        amax = zero;
        return;
    }
    //
    REAL base = Rlamch("B");
    REAL tmp = -0.5e+0 / log(base);
    //
    //     Find the minimum and maximum diagonal elements.
    //
    s[1 - 1] = a[(1 - 1)];
    REAL smin = s[1 - 1];
    amax = s[1 - 1];
    INTEGER i = 0;
    for (i = 2; i <= n; i = i + 1) {
        s[i - 1] = a[(i - 1) + (i - 1) * lda];
        smin = min(smin, s[i - 1]);
        amax = max(amax, s[i - 1]);
    }
    //
    if (smin <= zero) {
        //
        //        Find the first non-positive diagonal element and return.
        //
        for (i = 1; i <= n; i = i + 1) {
            if (s[i - 1] <= zero) {
                info = i;
                return;
            }
        }
    } else {
        //
        //        Set the scale factors to the reciprocals
        //        of the diagonal elements.
        //
        for (i = 1; i <= n; i = i + 1) {
            s[i - 1] = pow(base, int(tmp * log(s[i - 1])));
        }
        //
        //        Compute SCOND = min(S(I)) / max(S(I)).
        //
        scond = sqrt(smin) / sqrt(amax);
    }
    //
    //     End of Cpoequb
    //
}
