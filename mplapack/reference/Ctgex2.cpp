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

void Ctgex2(bool const wantq, bool const wantz, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *q, INTEGER const ldq, COMPLEX *z, INTEGER const ldz, INTEGER const j1, INTEGER &info) {
    const INTEGER ldst = 2;
    INTEGER m = 0;
    bool weak = false;
    bool strong = false;
    arr_2d<ldst, ldst, COMPLEX> s(fill0);
    arr_2d<ldst, ldst, COMPLEX> t(fill0);
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const COMPLEX czero = (0.0, 0.0);
    REAL scale = 0.0;
    const COMPLEX cone = (1.0, 0.0);
    REAL sum = 0.0;
    arr_1d<8, COMPLEX> work(fill0);
    REAL sa = 0.0;
    REAL sb = 0.0;
    const REAL twenty = 2.0e+1;
    REAL thresha = 0.0;
    REAL threshb = 0.0;
    COMPLEX f = 0.0;
    COMPLEX g = 0.0;
    REAL cz = 0.0;
    COMPLEX sz = 0.0;
    COMPLEX cdum = 0.0;
    REAL cq = 0.0;
    COMPLEX sq = 0.0;
    const bool wands = true;
    INTEGER i = 0;
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
    info = 0;
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    //
    m = ldst;
    weak = false;
    strong = false;
    //
    //     Make a local copy of selected block in (A, B)
    //
    Clacpy("Full", m, m, &a[(j1 - 1) + (j1 - 1) * lda], lda, s, ldst);
    Clacpy("Full", m, m, &b[(j1 - 1) + (j1 - 1) * ldb], ldb, t, ldst);
    //
    //     Compute the threshold for testing the acceptance of swapping.
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    scale = czero.real();
    sum = cone.real();
    Clacpy("Full", m, m, s, ldst, work, m);
    Clacpy("Full", m, m, t, ldst, &work[(m * m + 1) - 1], m);
    Classq(m * m, work, 1, scale, sum);
    sa = scale * sqrt(sum);
    scale = czero.real();
    sum = cone.real();
    Classq(m * m, &work[(m * m + 1) - 1], 1, scale, sum);
    sb = scale * sqrt(sum);
    //
    //     THRES has been changed from
    //        THRESH = MAX( TEN*EPS*SA, SMLNUM )
    //     to
    //        THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
    //     on 04/01/10.
    //     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
    //     Jim Demmel and Guillaume Revy. See forum post 1783.
    //
    thresha = max(twenty * eps * sa, smlnum);
    threshb = max(twenty * eps * sb, smlnum);
    //
    //     Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
    //     using Givens rotations and perform the swap tentatively.
    //
    f = s[(2 - 1) + (2 - 1) * lds] * t[(1 - 1)] - t[(2 - 1) + (2 - 1) * ldt] * s[(1 - 1)];
    g = s[(2 - 1) + (2 - 1) * lds] * t[(2 - 1) * ldt] - t[(2 - 1) + (2 - 1) * ldt] * s[(2 - 1) * lds];
    sa = abs(s[(2 - 1) + (2 - 1) * lds]) * abs(t[(1 - 1)]);
    sb = abs(s[(1 - 1)]) * abs(t[(2 - 1) + (2 - 1) * ldt]);
    Clartg(g, f, cz, sz, cdum);
    sz = -sz;
    Crot(2, s[(1 - 1)], 1, s[(2 - 1) * lds], 1, cz, conj(sz));
    Crot(2, &t[(1 - 1)], 1, &t[(2 - 1) * ldt], 1, cz, conj(sz));
    if (sa >= sb) {
        Clartg(s[(1 - 1)], s[(2 - 1)], cq, sq, cdum);
    } else {
        Clartg(t[(1 - 1)], &t[(2 - 1)], cq, sq, cdum);
    }
    Crot(2, s[(1 - 1)], ldst, s[(2 - 1)], ldst, cq, sq);
    Crot(2, &t[(1 - 1)], ldst, &t[(2 - 1)], ldst, cq, sq);
    //
    //     Weak stability test: |S21| <= O(EPS F-norm((A)))
    //                          and  |T21| <= O(EPS F-norm((B)))
    //
    weak = abs(s[(2 - 1)]) <= thresha && abs(t[(2 - 1)]) <= threshb;
    if (!weak) {
        goto statement_20;
    }
    //
    if (wands) {
        //
        //        Strong stability test:
        //           F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
        //           and
        //           F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
        //
        Clacpy("Full", m, m, s, ldst, work, m);
        Clacpy("Full", m, m, t, ldst, &work[(m * m + 1) - 1], m);
        Crot(2, work, 1, &work[3 - 1], 1, cz, -conj(sz));
        Crot(2, &work[5 - 1], 1, &work[7 - 1], 1, cz, -conj(sz));
        Crot(2, work, 2, &work[2 - 1], 2, cq, -sq);
        Crot(2, &work[5 - 1], 2, &work[6 - 1], 2, cq, -sq);
        for (i = 1; i <= 2; i = i + 1) {
            work[i - 1] = work[i - 1] - a[((j1 + i - 1) - 1) + (j1 - 1) * lda];
            work[(i + 2) - 1] = work[(i + 2) - 1] - a[((j1 + i - 1) - 1) + ((j1 + 1) - 1) * lda];
            work[(i + 4) - 1] = work[(i + 4) - 1] - b[((j1 + i - 1) - 1) + (j1 - 1) * ldb];
            work[(i + 6) - 1] = work[(i + 6) - 1] - b[((j1 + i - 1) - 1) + ((j1 + 1) - 1) * ldb];
        }
        scale = czero.real();
        sum = cone.real();
        Classq(m * m, work, 1, scale, sum);
        sa = scale * sqrt(sum);
        scale = czero.real();
        sum = cone.real();
        Classq(m * m, &work[(m * m + 1) - 1], 1, scale, sum);
        sb = scale * sqrt(sum);
        strong = sa <= thresha && sb <= threshb;
        if (!strong) {
            goto statement_20;
        }
    }
    //
    //     If the swap is accepted ("weakly" and "strongly"), apply the
    //     equivalence transformations to the original matrix pair (A,B)
    //
    Crot(j1 + 1, &a[(j1 - 1) * lda], 1, &a[((j1 + 1) - 1) * lda], 1, cz, conj(sz));
    Crot(j1 + 1, &b[(j1 - 1) * ldb], 1, &b[((j1 + 1) - 1) * ldb], 1, cz, conj(sz));
    Crot(n - j1 + 1, &a[(j1 - 1) + (j1 - 1) * lda], lda, &a[((j1 + 1) - 1) + (j1 - 1) * lda], lda, cq, sq);
    Crot(n - j1 + 1, &b[(j1 - 1) + (j1 - 1) * ldb], ldb, &b[((j1 + 1) - 1) + (j1 - 1) * ldb], ldb, cq, sq);
    //
    //     Set  N1 by N2 (2,1) blocks to 0
    //
    a[((j1 + 1) - 1) + (j1 - 1) * lda] = czero;
    b[((j1 + 1) - 1) + (j1 - 1) * ldb] = czero;
    //
    //     Accumulate transformations into Q and Z if requested.
    //
    if (wantz) {
        Crot(n, &z[(j1 - 1) * ldz], 1, &z[((j1 + 1) - 1) * ldz], 1, cz, conj(sz));
    }
    if (wantq) {
        Crot(n, q[(j1 - 1) * ldq], 1, q[((j1 + 1) - 1) * ldq], 1, cq, conj(sq));
    }
    //
    //     Exit with INFO = 0 if swap was successfully performed.
    //
    return;
//
//     Exit with INFO = 1 if swap was rejected.
//
statement_20:
    info = 1;
    //
    //     End of Ctgex2
    //
}
