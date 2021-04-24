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

void Rlaror(const char *side, const char *init, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *iseed, REAL *x, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if (n == 0 || m == 0) {
        return;
    }
    //
    INTEGER itype = 0;
    if (Mlsame(side, "L")) {
        itype = 1;
    } else if (Mlsame(side, "R")) {
        itype = 2;
    } else if (Mlsame(side, "C") || Mlsame(side, "T")) {
        itype = 3;
    }
    //
    //     Check for argument errors.
    //
    if (itype == 0) {
        info = -1;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0 || (itype == 3 && n != m)) {
        info = -4;
    } else if (lda < m) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rlaror", -info);
        return;
    }
    //
    INTEGER nxfrm = 0;
    if (itype == 1) {
        nxfrm = m;
    } else {
        nxfrm = n;
    }
    //
    //     Initialize A to the identity matrix if desired
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (Mlsame(init, "I")) {
        dlaset("Full", m, n, zero, one, a, lda);
    }
    //
    //     If no rotation possible, multiply by random +/-1
    //
    //     Compute rotation by computing Householder transformations
    //     H(2), H(3), ..., H(nhouse)
    //
    INTEGER j = 0;
    for (j = 1; j <= nxfrm; j = j + 1) {
        x[j - 1] = zero;
    }
    //
    INTEGER ixfrm = 0;
    INTEGER kbeg = 0;
    REAL xnorm = 0.0;
    REAL xnorms = 0.0;
    REAL factor = 0.0;
    const REAL toosml = 1.0e-20;
    for (ixfrm = 2; ixfrm <= nxfrm; ixfrm = ixfrm + 1) {
        kbeg = nxfrm - ixfrm + 1;
        //
        //        Generate independent normal( 0, 1 ) random numbers
        //
        for (j = kbeg; j <= nxfrm; j = j + 1) {
            x[j - 1] = Rlarnd[(3 - 1) + (iseed - 1) * ldRlarnd];
        }
        //
        //        Generate a Householder transformation from the random vector X
        //
        xnorm = Rnrm2(ixfrm, &x[kbeg - 1], 1);
        xnorms = sign(xnorm, &x[kbeg - 1]);
        x[(kbeg + nxfrm) - 1] = sign(one, -x[kbeg - 1]);
        factor = xnorms * (xnorms + x[kbeg - 1]);
        if (abs(factor) < toosml) {
            info = 1;
            Mxerbla("Rlaror", info);
            return;
        } else {
            factor = one / factor;
        }
        x[kbeg - 1] += xnorms;
        //
        //        Apply Householder transformation to A
        //
        if (itype == 1 || itype == 3) {
            //
            //           Apply H(k) from the left.
            //
            Rgemv("T", ixfrm, n, one, &a[(kbeg - 1)], lda, &x[kbeg - 1], 1, zero, &x[(2 * nxfrm + 1) - 1], 1);
            Rger(ixfrm, n, -factor, &x[kbeg - 1], 1, &x[(2 * nxfrm + 1) - 1], 1, &a[(kbeg - 1)], lda);
            //
        }
        //
        if (itype == 2 || itype == 3) {
            //
            //           Apply H(k) from the right.
            //
            Rgemv("N", m, ixfrm, one, &a[(kbeg - 1) * lda], lda, &x[kbeg - 1], 1, zero, &x[(2 * nxfrm + 1) - 1], 1);
            Rger(m, ixfrm, -factor, &x[(2 * nxfrm + 1) - 1], 1, &x[kbeg - 1], 1, &a[(kbeg - 1) * lda], lda);
            //
        }
    }
    //
    x[(2 * nxfrm) - 1] = sign(one, Rlarnd[(3 - 1) + (iseed - 1) * ldRlarnd]);
    //
    //     Scale the matrix A by D.
    //
    INTEGER irow = 0;
    if (itype == 1 || itype == 3) {
        for (irow = 1; irow <= m; irow = irow + 1) {
            Rscal(n, &x[(nxfrm + irow) - 1], &a[(irow - 1)], lda);
        }
    }
    //
    INTEGER jcol = 0;
    if (itype == 2 || itype == 3) {
        for (jcol = 1; jcol <= n; jcol = jcol + 1) {
            Rscal(m, &x[(nxfrm + jcol) - 1], &a[(jcol - 1) * lda], 1);
        }
    }
    //
    //     End of Rlaror
    //
}
