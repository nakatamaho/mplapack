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

void Rtgex2(bool const &wantq, bool const &wantz, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *q, INTEGER const &ldq, REAL *z, INTEGER const &ldz, INTEGER const &j1, INTEGER const &n1, INTEGER const &n2, REAL *work, INTEGER const &lwork, INTEGER &info) {
    INTEGER m = 0;
    bool weak = false;
    bool strong = false;
    const INTEGER ldst = 4;
    const REAL zero = 0.0;
    arr_2d<ldst, ldst, REAL> li(fill0);
    arr_2d<ldst, ldst, REAL> ir(fill0);
    arr_2d<ldst, ldst, REAL> s(fill0);
    arr_2d<ldst, ldst, REAL> t(fill0);
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL Rscale = 0.0;
    const REAL one = 1.0;
    REAL dsum = 0.0;
    REAL dnorma = 0.0;
    REAL dnormb = 0.0;
    const REAL twenty = 2.0e+01;
    REAL thresha = 0.0;
    REAL threshb = 0.0;
    REAL f = 0.0;
    REAL g = 0.0;
    REAL sa = 0.0;
    REAL sb = 0.0;
    REAL ddum = 0.0;
    const bool wands = true;
    REAL scale = 0.0;
    arr_1d<ldst, INTEGER> iwork(fill0);
    INTEGER idum = 0;
    INTEGER linfo = 0;
    INTEGER i = 0;
    arr_1d<ldst, REAL> taul(fill0);
    arr_1d<ldst, REAL> taur(fill0);
    arr_2d<ldst, ldst, REAL> scpy(fill0);
    arr_2d<ldst, ldst, REAL> tcpy(fill0);
    arr_2d<ldst, ldst, REAL> ircop(fill0);
    arr_2d<ldst, ldst, REAL> licop(fill0);
    REAL brqa21 = 0.0;
    REAL bqra21 = 0.0;
    arr_1d<2, REAL> ar(fill0);
    arr_1d<2, REAL> ai(fill0);
    arr_1d<2, REAL> be(fill0);
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
    //  Replaced various illegal calls to Rcopy by calls to Rlaset, or by DO
    //  loops. Sven Hammarling, 1/5/02.
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
    if (n <= 1 || n1 <= 0 || n2 <= 0) {
        return;
    }
    if (n1 > n || (j1 + n1) > n) {
        return;
    }
    m = n1 + n2;
    if (lwork < max((INTEGER)1, n * m, m * m * 2)) {
        info = -16;
        work[1 - 1] = max((INTEGER)1, n * m, m * m * 2);
        return;
    }
    //
    weak = false;
    strong = false;
    //
    //     Make a local copy of selected block
    //
    Rlaset("Full", ldst, ldst, zero, zero, li, ldst);
    Rlaset("Full", ldst, ldst, zero, zero, ir, ldst);
    Rlacpy("Full", m, m, a[(j1 - 1) + (j1 - 1) * lda], lda, s, ldst);
    Rlacpy("Full", m, m, b[(j1 - 1) + (j1 - 1) * ldb], ldb, t, ldst);
    //
    //     Compute threshold for testing acceptance of swapping.
    //
    eps = dlamch("P");
    smlnum = dlamch("S") / eps;
    Rscale = zero;
    dsum = one;
    Rlacpy("Full", m, m, s, ldst, work, m);
    Rlassq(m * m, work, 1, Rscale, dsum);
    dnorma = Rscale * sqrt(dsum);
    Rscale = zero;
    dsum = one;
    Rlacpy("Full", m, m, t, ldst, work, m);
    Rlassq(m * m, work, 1, Rscale, dsum);
    dnormb = Rscale * sqrt(dsum);
    //
    //     THRES has been changed from
    //        THRESH = MAX( TEN*EPS*SA, SMLNUM )
    //     to
    //        THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
    //     on 04/01/10.
    //     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
    //     Jim Demmel and Guillaume Revy. See forum post 1783.
    //
    thresha = max(twenty * eps * dnorma, smlnum);
    threshb = max(twenty * eps * dnormb, smlnum);
    //
    if (m == 2) {
        //
        //        CASE 1: Swap 1-by-1 and 1-by-1 blocks.
        //
        //        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
        //        using Givens rotations and perform the swap tentatively.
        //
        f = s[(2 - 1) + (2 - 1) * lds] * t[(1 - 1)] - t[(2 - 1) + (2 - 1) * ldt] * s[(1 - 1)];
        g = s[(2 - 1) + (2 - 1) * lds] * t[(2 - 1) * ldt] - t[(2 - 1) + (2 - 1) * ldt] * s[(2 - 1) * lds];
        sa = abs(s[(2 - 1) + (2 - 1) * lds]) * abs(t[(1 - 1)]);
        sb = abs(s[(1 - 1)]) * abs(t[(2 - 1) + (2 - 1) * ldt]);
        Rlartg(f, g, ir[(2 - 1) * ldir], ir[(1 - 1)], ddum);
        ir[(2 - 1)] = -ir[(2 - 1) * ldir];
        ir[(2 - 1) + (2 - 1) * ldir] = ir[(1 - 1)];
        Rrot(2, s[(1 - 1)], 1, s[(2 - 1) * lds], 1, ir[(1 - 1)], ir[(2 - 1)]);
        Rrot(2, t[(1 - 1)], 1, t[(2 - 1) * ldt], 1, ir[(1 - 1)], ir[(2 - 1)]);
        if (sa >= sb) {
            Rlartg(s[(1 - 1)], s[(2 - 1)], li[(1 - 1)], li[(2 - 1)], ddum);
        } else {
            Rlartg(t[(1 - 1)], t[(2 - 1)], li[(1 - 1)], li[(2 - 1)], ddum);
        }
        Rrot(2, s[(1 - 1)], ldst, s[(2 - 1)], ldst, li[(1 - 1)], li[(2 - 1)]);
        Rrot(2, t[(1 - 1)], ldst, t[(2 - 1)], ldst, li[(1 - 1)], li[(2 - 1)]);
        li[(2 - 1) + (2 - 1) * ldli] = li[(1 - 1)];
        li[(2 - 1) * ldli] = -li[(2 - 1)];
        //
        //        Weak stability test: |S21| <= O(EPS F-norm((A)))
        //                           and  |T21| <= O(EPS F-norm((B)))
        //
        weak = abs(s[(2 - 1)]) <= thresha && abs(t[(2 - 1)]) <= threshb;
        if (!weak) {
            goto statement_70;
        }
        //
        if (wands) {
            //
            //           Strong stability test:
            //               F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
            //               and
            //               F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
            //
            Rlacpy("Full", m, m, a[(j1 - 1) + (j1 - 1) * lda], lda, work[(m * m + 1) - 1], m);
            Rgemm("N", "N", m, m, m, one, li, ldst, s, ldst, zero, work, m);
            Rgemm("N", "T", m, m, m, -one, work, m, ir, ldst, one, work[(m * m + 1) - 1], m);
            Rscale = zero;
            dsum = one;
            Rlassq(m * m, work[(m * m + 1) - 1], 1, Rscale, dsum);
            sa = Rscale * sqrt(dsum);
            //
            Rlacpy("Full", m, m, b[(j1 - 1) + (j1 - 1) * ldb], ldb, work[(m * m + 1) - 1], m);
            Rgemm("N", "N", m, m, m, one, li, ldst, t, ldst, zero, work, m);
            Rgemm("N", "T", m, m, m, -one, work, m, ir, ldst, one, work[(m * m + 1) - 1], m);
            Rscale = zero;
            dsum = one;
            Rlassq(m * m, work[(m * m + 1) - 1], 1, Rscale, dsum);
            sb = Rscale * sqrt(dsum);
            strong = sa <= thresha && sb <= threshb;
            if (!strong) {
                goto statement_70;
            }
        }
        //
        //        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
        //               (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
        //
        Rrot(j1 + 1, a[(j1 - 1) * lda], 1, a[((j1 + 1) - 1) * lda], 1, ir[(1 - 1)], ir[(2 - 1)]);
        Rrot(j1 + 1, b[(j1 - 1) * ldb], 1, b[((j1 + 1) - 1) * ldb], 1, ir[(1 - 1)], ir[(2 - 1)]);
        Rrot(n - j1 + 1, a[(j1 - 1) + (j1 - 1) * lda], lda, a[((j1 + 1) - 1) + (j1 - 1) * lda], lda, li[(1 - 1)], li[(2 - 1)]);
        Rrot(n - j1 + 1, b[(j1 - 1) + (j1 - 1) * ldb], ldb, b[((j1 + 1) - 1) + (j1 - 1) * ldb], ldb, li[(1 - 1)], li[(2 - 1)]);
        //
        //        Set  N1-by-N2 (2,1) - blocks to ZERO.
        //
        a[((j1 + 1) - 1) + (j1 - 1) * lda] = zero;
        b[((j1 + 1) - 1) + (j1 - 1) * ldb] = zero;
        //
        //        Accumulate transformations INTEGERo Q and Z if requested.
        //
        if (wantz) {
            Rrot(n, z[(j1 - 1) * ldz], 1, z[((j1 + 1) - 1) * ldz], 1, ir[(1 - 1)], ir[(2 - 1)]);
        }
        if (wantq) {
            Rrot(n, q[(j1 - 1) * ldq], 1, q[((j1 + 1) - 1) * ldq], 1, li[(1 - 1)], li[(2 - 1)]);
        }
        //
        //        Exit with INFO = 0 if swap was successfully performed.
        //
        return;
        //
    } else {
        //
        //        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
        //                and 2-by-2 blocks.
        //
        //        Solve the generalized Sylvester equation
        //                 S11 * R - L * S22 = SCALE * S12
        //                 T11 * R - L * T22 = SCALE * T12
        //        for R and L. Solutions in LI and IR.
        //
        Rlacpy("Full", n1, n2, t[((n1 + 1) - 1) * ldt], ldst, li, ldst);
        Rlacpy("Full", n1, n2, s[((n1 + 1) - 1) * lds], ldst, ir[((n2 + 1) - 1) + ((n1 + 1) - 1) * ldir], ldst);
        Rtgsy2("N", 0, n1, n2, s, ldst, s[((n1 + 1) - 1) + ((n1 + 1) - 1) * lds], ldst, ir[((n2 + 1) - 1) + ((n1 + 1) - 1) * ldir], ldst, t, ldst, t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldst, li, ldst, scale, dsum, Rscale, iwork, idum, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        //
        //        Compute orthogonal matrix QL:
        //
        //                    QL**T * LI = [ TL ]
        //                                 [ 0  ]
        //        where
        //                    LI =  [      -L              ]
        //                          [ SCALE * identity(N2) ]
        //
        for (i = 1; i <= n2; i = i + 1) {
            Rscal(n1, -one, li[(i - 1) * ldli], 1);
            li[((n1 + i) - 1) + (i - 1) * ldli] = scale;
        }
        Rgeqr2(m, n2, li, ldst, taul, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        Rorg2r(m, m, n2, li, ldst, taul, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        //
        //        Compute orthogonal matrix RQ:
        //
        //                    IR * RQ**T =   [ 0  TR],
        //
        //         where IR = [ SCALE * identity(N1), R ]
        //
        for (i = 1; i <= n1; i = i + 1) {
            ir[((n2 + i) - 1) + (i - 1) * ldir] = scale;
        }
        Rgerq2(n1, m, ir[((n2 + 1) - 1)], ldst, taur, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        Rorgr2(m, m, n1, ir, ldst, taur, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        //
        //        Perform the swapping tentatively:
        //
        Rgemm("T", "N", m, m, m, one, li, ldst, s, ldst, zero, work, m);
        Rgemm("N", "T", m, m, m, one, work, m, ir, ldst, zero, s, ldst);
        Rgemm("T", "N", m, m, m, one, li, ldst, t, ldst, zero, work, m);
        Rgemm("N", "T", m, m, m, one, work, m, ir, ldst, zero, t, ldst);
        Rlacpy("F", m, m, s, ldst, scpy, ldst);
        Rlacpy("F", m, m, t, ldst, tcpy, ldst);
        Rlacpy("F", m, m, ir, ldst, ircop, ldst);
        Rlacpy("F", m, m, li, ldst, licop, ldst);
        //
        //        Triangularize the B-part by an RQ factorization.
        //        Apply transformation (from left) to A-part, giving S.
        //
        Rgerq2(m, m, t, ldst, taur, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        Rormr2("R", "T", m, m, m, t, ldst, taur, s, ldst, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        Rormr2("L", "N", m, m, m, t, ldst, taur, ir, ldst, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        //
        //        Compute F-norm(S21) in BRQA21. (T21 is 0.)
        //
        Rscale = zero;
        dsum = one;
        for (i = 1; i <= n2; i = i + 1) {
            Rlassq(n1, s[((n2 + 1) - 1) + (i - 1) * lds], 1, Rscale, dsum);
        }
        brqa21 = Rscale * sqrt(dsum);
        //
        //        Triangularize the B-part by a QR factorization.
        //        Apply transformation (from right) to A-part, giving S.
        //
        Rgeqr2(m, m, tcpy, ldst, taul, work, linfo);
        if (linfo != 0) {
            goto statement_70;
        }
        Rorm2r("L", "T", m, m, m, tcpy, ldst, taul, scpy, ldst, work, info);
        Rorm2r("R", "N", m, m, m, tcpy, ldst, taul, licop, ldst, work, info);
        if (linfo != 0) {
            goto statement_70;
        }
        //
        //        Compute F-norm(S21) in BQRA21. (T21 is 0.)
        //
        Rscale = zero;
        dsum = one;
        for (i = 1; i <= n2; i = i + 1) {
            Rlassq(n1, scpy[((n2 + 1) - 1) + (i - 1) * ldscpy], 1, Rscale, dsum);
        }
        bqra21 = Rscale * sqrt(dsum);
        //
        //        Decide which method to use.
        //          Weak stability test:
        //             F-norm(S21) <= O(EPS * F-norm((S)))
        //
        if (bqra21 <= brqa21 && bqra21 <= thresha) {
            Rlacpy("F", m, m, scpy, ldst, s, ldst);
            Rlacpy("F", m, m, tcpy, ldst, t, ldst);
            Rlacpy("F", m, m, ircop, ldst, ir, ldst);
            Rlacpy("F", m, m, licop, ldst, li, ldst);
        } else if (brqa21 >= thresha) {
            goto statement_70;
        }
        //
        //        Set lower triangle of B-part to zero
        //
        Rlaset("Lower", m - 1, m - 1, zero, zero, t[(2 - 1)], ldst);
        //
        if (wands) {
            //
            //           Strong stability test:
            //               F-norm((A-QL**H*S*QR)) <= O(EPS*F-norm((A)))
            //               and
            //               F-norm((B-QL**H*T*QR)) <= O(EPS*F-norm((B)))
            //
            Rlacpy("Full", m, m, a[(j1 - 1) + (j1 - 1) * lda], lda, work[(m * m + 1) - 1], m);
            Rgemm("N", "N", m, m, m, one, li, ldst, s, ldst, zero, work, m);
            Rgemm("N", "N", m, m, m, -one, work, m, ir, ldst, one, work[(m * m + 1) - 1], m);
            Rscale = zero;
            dsum = one;
            Rlassq(m * m, work[(m * m + 1) - 1], 1, Rscale, dsum);
            sa = Rscale * sqrt(dsum);
            //
            Rlacpy("Full", m, m, b[(j1 - 1) + (j1 - 1) * ldb], ldb, work[(m * m + 1) - 1], m);
            Rgemm("N", "N", m, m, m, one, li, ldst, t, ldst, zero, work, m);
            Rgemm("N", "N", m, m, m, -one, work, m, ir, ldst, one, work[(m * m + 1) - 1], m);
            Rscale = zero;
            dsum = one;
            Rlassq(m * m, work[(m * m + 1) - 1], 1, Rscale, dsum);
            sb = Rscale * sqrt(dsum);
            strong = sa <= thresha && sb <= threshb;
            if (!strong) {
                goto statement_70;
            }
            //
        }
        //
        //        If the swap is accepted ("weakly" and "strongly"), apply the
        //        transformations and set N1-by-N2 (2,1)-block to zero.
        //
        Rlaset("Full", n1, n2, zero, zero, s[((n2 + 1) - 1)], ldst);
        //
        //        copy back M-by-M diagonal block starting at index J1 of (A, B)
        //
        Rlacpy("F", m, m, s, ldst, a[(j1 - 1) + (j1 - 1) * lda], lda);
        Rlacpy("F", m, m, t, ldst, b[(j1 - 1) + (j1 - 1) * ldb], ldb);
        Rlaset("Full", ldst, ldst, zero, zero, t, ldst);
        //
        //        Standardize existing 2-by-2 blocks.
        //
        Rlaset("Full", m, m, zero, zero, work, m);
        work[1 - 1] = one;
        t[(1 - 1)] = one;
        idum = lwork - m * m - 2;
        if (n2 > 1) {
            Rlagv2(a[(j1 - 1) + (j1 - 1) * lda], lda, b[(j1 - 1) + (j1 - 1) * ldb], ldb, ar, ai, be, work[1 - 1], work[2 - 1], t[(1 - 1)], t[(2 - 1)]);
            work[(m + 1) - 1] = -work[2 - 1];
            work[(m + 2) - 1] = work[1 - 1];
            t[(n2 - 1) + (n2 - 1) * ldt] = t[(1 - 1)];
            t[(2 - 1) * ldt] = -t[(2 - 1)];
        }
        work[(m * m) - 1] = one;
        t[(m - 1) + (m - 1) * ldt] = one;
        //
        if (n1 > 1) {
            Rlagv2(a[((j1 + n2) - 1) + ((j1 + n2) - 1) * lda], lda, b[((j1 + n2) - 1) + ((j1 + n2) - 1) * ldb], ldb, taur, taul, work[(m * m + 1) - 1], work[(n2 * m + n2 + 1) - 1], work[(n2 * m + n2 + 2) - 1], t[((n2 + 1) - 1) + ((n2 + 1) - 1) * ldt], t[(m - 1) + ((m - 1) - 1) * ldt]);
            work[(m * m) - 1] = work[(n2 * m + n2 + 1) - 1];
            work[(m * m - 1) - 1] = -work[(n2 * m + n2 + 2) - 1];
            t[(m - 1) + (m - 1) * ldt] = t[((n2 + 1) - 1) + ((n2 + 1) - 1) * ldt];
            t[((m - 1) - 1) + (m - 1) * ldt] = -t[(m - 1) + ((m - 1) - 1) * ldt];
        }
        Rgemm("T", "N", n2, n1, n2, one, work, m, a[(j1 - 1) + ((j1 + n2) - 1) * lda], lda, zero, work[(m * m + 1) - 1], n2);
        Rlacpy("Full", n2, n1, work[(m * m + 1) - 1], n2, a[(j1 - 1) + ((j1 + n2) - 1) * lda], lda);
        Rgemm("T", "N", n2, n1, n2, one, work, m, b[(j1 - 1) + ((j1 + n2) - 1) * ldb], ldb, zero, work[(m * m + 1) - 1], n2);
        Rlacpy("Full", n2, n1, work[(m * m + 1) - 1], n2, b[(j1 - 1) + ((j1 + n2) - 1) * ldb], ldb);
        Rgemm("N", "N", m, m, m, one, li, ldst, work, m, zero, work[(m * m + 1) - 1], m);
        Rlacpy("Full", m, m, work[(m * m + 1) - 1], m, li, ldst);
        Rgemm("N", "N", n2, n1, n1, one, a[(j1 - 1) + ((j1 + n2) - 1) * lda], lda, t[((n2 + 1) - 1) + ((n2 + 1) - 1) * ldt], ldst, zero, work, n2);
        Rlacpy("Full", n2, n1, work, n2, a[(j1 - 1) + ((j1 + n2) - 1) * lda], lda);
        Rgemm("N", "N", n2, n1, n1, one, b[(j1 - 1) + ((j1 + n2) - 1) * ldb], ldb, t[((n2 + 1) - 1) + ((n2 + 1) - 1) * ldt], ldst, zero, work, n2);
        Rlacpy("Full", n2, n1, work, n2, b[(j1 - 1) + ((j1 + n2) - 1) * ldb], ldb);
        Rgemm("T", "N", m, m, m, one, ir, ldst, t, ldst, zero, work, m);
        Rlacpy("Full", m, m, work, m, ir, ldst);
        //
        //        Accumulate transformations INTEGERo Q and Z if requested.
        //
        if (wantq) {
            Rgemm("N", "N", n, m, m, one, q[(j1 - 1) * ldq], ldq, li, ldst, zero, work, n);
            Rlacpy("Full", n, m, work, n, q[(j1 - 1) * ldq], ldq);
            //
        }
        //
        if (wantz) {
            Rgemm("N", "N", n, m, m, one, z[(j1 - 1) * ldz], ldz, ir, ldst, zero, work, n);
            Rlacpy("Full", n, m, work, n, z[(j1 - 1) * ldz], ldz);
            //
        }
        //
        //        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
        //                (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
        //
        i = j1 + m;
        if (i <= n) {
            Rgemm("T", "N", m, n - i + 1, m, one, li, ldst, a[(j1 - 1) + (i - 1) * lda], lda, zero, work, m);
            Rlacpy("Full", m, n - i + 1, work, m, a[(j1 - 1) + (i - 1) * lda], lda);
            Rgemm("T", "N", m, n - i + 1, m, one, li, ldst, b[(j1 - 1) + (i - 1) * ldb], ldb, zero, work, m);
            Rlacpy("Full", m, n - i + 1, work, m, b[(j1 - 1) + (i - 1) * ldb], ldb);
        }
        i = j1 - 1;
        if (i > 0) {
            Rgemm("N", "N", i, m, m, one, a[(j1 - 1) * lda], lda, ir, ldst, zero, work, i);
            Rlacpy("Full", i, m, work, i, a[(j1 - 1) * lda], lda);
            Rgemm("N", "N", i, m, m, one, b[(j1 - 1) * ldb], ldb, ir, ldst, zero, work, i);
            Rlacpy("Full", i, m, work, i, b[(j1 - 1) * ldb], ldb);
        }
        //
        //        Exit with INFO = 0 if swap was successfully performed.
        //
        return;
        //
    }
//
//     Exit with INFO = 1 if swap was rejected.
//
statement_70:
    //
    info = 1;
    //
    //     End of Rtgex2
    //
}
