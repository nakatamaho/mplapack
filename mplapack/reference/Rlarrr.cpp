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

void Rlarrr(INTEGER const n, REAL *d, REAL *e, INTEGER &info) {
    REAL safmin = 0.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL rmin = 0.0;
    bool yesrel = false;
    const REAL zero = 0.0;
    REAL offdig = 0.0;
    REAL tmp = 0.0;
    INTEGER i = 0;
    REAL tmp2 = 0.0;
    REAL offdig2 = 0.0;
    const REAL relcond = 0.999e0;
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
    //
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        info = 0;
        return;
    }
    //
    //     As a default, do NOT go for relative-accuracy preserving computations.
    info = 1;
    //
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
    smlnum = safmin / eps;
    rmin = sqrt(smlnum);
    //
    //     Tests for relative accuracy
    //
    //     Test for scaled diagonal dominance
    //     Scale the diagonal entries to one and check whether the sum of the
    //     off-diagonals is less than one
    //
    //     The sdd relative error bounds have a 1/(1- 2*x) factor in them,
    //     x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative
    //     accuracy is promised.  In the notation of the code fragment below,
    //     1/(1 - (OFFDIG + OFFDIG2)) is the condition number.
    //     We don't think it is worth going into "sdd mode" unless the relative
    //     condition number is reasonable, not 1/macheps.
    //     The threshold should be compatible with other thresholds used in the
    //     code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds
    //     to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000
    //     instead of the current OFFDIG + OFFDIG2 < 1
    //
    yesrel = true;
    offdig = zero;
    tmp = sqrt(abs(d[1 - 1]));
    if (tmp < rmin) {
        yesrel = false;
    }
    if (!yesrel) {
        goto statement_11;
    }
    for (i = 2; i <= n; i = i + 1) {
        tmp2 = sqrt(abs(d[i - 1]));
        if (tmp2 < rmin) {
            yesrel = false;
        }
        if (!yesrel) {
            goto statement_11;
        }
        offdig2 = abs(e[(i - 1) - 1]) / (tmp * tmp2);
        if (offdig + offdig2 >= relcond) {
            yesrel = false;
        }
        if (!yesrel) {
            goto statement_11;
        }
        tmp = tmp2;
        offdig = offdig2;
    }
statement_11:
    //
    if (yesrel) {
        info = 0;
        return;
    } else {
    }
    //
    //     *** MORE TO BE IMPLEMENTED ***
    //
    //     Test if the lower bidiagonal matrix L from T = L D L^T
    //     (zero shift facto) is well conditioned
    //
    //     Test if the upper bidiagonal matrix U from T = U D U^T
    //     (zero shift facto) is well conditioned.
    //     In this case, the matrix needs to be flipped and, at the end
    //     of the eigenvector computation, the flip needs to be applied
    //     to the computed eigenvectors (and the support)
    //
    //     END OF Rlarrr
    //
}
