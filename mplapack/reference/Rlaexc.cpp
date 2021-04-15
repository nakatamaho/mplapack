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

void Rlaexc(bool const wantq, INTEGER const n, REAL *t, INTEGER const ldt, REAL *q, INTEGER const ldq, INTEGER const j1, INTEGER const n1, INTEGER const n2, REAL *work, INTEGER &info) {
    INTEGER j2 = 0;
    INTEGER j3 = 0;
    INTEGER j4 = 0;
    REAL t11 = 0.0;
    REAL t22 = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL temp = 0.0;
    INTEGER nd = 0;
    const INTEGER ldd = 4;
    REAL d[ldd * 4];
    REAL dnorm = 0.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL ten = 1.0e+1;
    REAL thresh = 0.0;
    REAL scale = 0.0;
    const INTEGER ldx = 2;
    REAL x[ldx * 2];
    REAL xnorm = 0.0;
    INTEGER ierr = 0;
    INTEGER k = 0;
    REAL u[3];
    REAL tau = 0.0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    REAL t33 = 0.0;
    REAL u1[3];
    REAL tau1 = 0.0;
    REAL u2[3];
    REAL tau2 = 0.0;
    REAL wr1 = 0.0;
    REAL wi1 = 0.0;
    REAL wr2 = 0.0;
    REAL wi2 = 0.0;
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
    if (n == 0 || n1 == 0 || n2 == 0) {
        return;
    }
    if (j1 + n1 > n) {
        return;
    }
    //
    j2 = j1 + 1;
    j3 = j1 + 2;
    j4 = j1 + 3;
    //
    if (n1 == 1 && n2 == 1) {
        //
        //        Swap two 1-by-1 blocks.
        //
        t11 = t[(j1 - 1) + (j1 - 1) * ldt];
        t22 = t[(j2 - 1) + (j2 - 1) * ldt];
        //
        //        Determine the transformation to perform the interchange.
        //
        Rlartg(t[(j1 - 1) + (j2 - 1) * ldt], t22 - t11, cs, sn, temp);
        //
        //        Apply transformation to the matrix T.
        //
        if (j3 <= n) {
            Rrot(n - j1 - 1, &t[(j1 - 1) + (j3 - 1) * ldt], ldt, &t[(j2 - 1) + (j3 - 1) * ldt], ldt, cs, sn);
        }
        Rrot(j1 - 1, &t[(j1 - 1) * ldt], 1, &t[(j2 - 1) * ldt], 1, cs, sn);
        //
        t[(j1 - 1) + (j1 - 1) * ldt] = t22;
        t[(j2 - 1) + (j2 - 1) * ldt] = t11;
        //
        if (wantq) {
            //
            //           Accumulate transformation in the matrix Q.
            //
            Rrot(n, &q[(j1 - 1) * ldq], 1, &q[(j2 - 1) * ldq], 1, cs, sn);
        }
        //
    } else {
        //
        //        Swapping involves at least one 2-by-2 block.
        //
        //        Copy the diagonal block of order N1+N2 to the local array D
        //        and compute its norm.
        //
        nd = n1 + n2;
        Rlacpy("Full", nd, nd, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, d, ldd);
        dnorm = Rlange("Max", nd, nd, d, ldd, work);
        //
        //        Compute machine-dependent threshold for test for accepting
        //        swap.
        //
        eps = Rlamch("P");
        smlnum = Rlamch("S") / eps;
        thresh = max(ten * eps * dnorm, smlnum);
        //
        //        Solve T11*X - X*T22 = scale*T12 for X.
        //
        Rlasy2(false, false, -1, n1, n2, d, ldd, &d[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldd], ldd, &d[((n1 + 1) - 1) * ldd], ldd, scale, x, ldx, xnorm, ierr);
        //
        //        Swap the adjacent diagonal blocks.
        //
        k = n1 + n1 + n2 - 3;
        switch (k) {
        case 1:
            goto statement_10;
        case 2:
            goto statement_20;
        case 3:
            goto statement_30;
        default:
            break;
        }
    //
    statement_10:
        //
        //        N1 = 1, N2 = 2: generate elementary reflector H so that:
        //
        //        ( scale, X11, X12 ) H = ( 0, 0, * )
        //
        u[1 - 1] = scale;
        u[2 - 1] = x[(1 - 1)];
        u[3 - 1] = x[(2 - 1) * ldx];
        Rlarfg(3, u[3 - 1], u, 1, tau);
        u[3 - 1] = one;
        t11 = t[(j1 - 1) + (j1 - 1) * ldt];
        //
        //        Perform swap provisionally on diagonal block in D.
        //
        Rlarfx("L", 3, 3, u, tau, d, ldd, work);
        Rlarfx("R", 3, 3, u, tau, d, ldd, work);
        //
        //        Test whether to reject swap.
        //
        if (max({abs(d[(3 - 1)]), abs(d[(3 - 1) + (2 - 1) * ldd]), abs(d[(3 - 1) + (3 - 1) * ldd] - t11)}) > thresh) {
            goto statement_50;
        }
        //
        //        Accept swap: apply transformation to the entire matrix T.
        //
        Rlarfx("L", 3, n - j1 + 1, u, tau, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, work);
        Rlarfx("R", j2, 3, u, tau, &t[(j1 - 1) * ldt], ldt, work);
        //
        t[(j3 - 1) + (j1 - 1) * ldt] = zero;
        t[(j3 - 1) + (j2 - 1) * ldt] = zero;
        t[(j3 - 1) + (j3 - 1) * ldt] = t11;
        //
        if (wantq) {
            //
            //           Accumulate transformation in the matrix Q.
            //
            Rlarfx("R", n, 3, u, tau, &q[(j1 - 1) * ldq], ldq, work);
        }
        goto statement_40;
    //
    statement_20:
        //
        //        N1 = 2, N2 = 1: generate elementary reflector H so that:
        //
        //        H (  -X11 ) = ( * )
        //          (  -X21 ) = ( 0 )
        //          ( scale ) = ( 0 )
        //
        u[1 - 1] = -x[(1 - 1)];
        u[2 - 1] = -x[(2 - 1)];
        u[3 - 1] = scale;
        Rlarfg(3, u[1 - 1], &u[2 - 1], 1, tau);
        u[1 - 1] = one;
        t33 = t[(j3 - 1) + (j3 - 1) * ldt];
        //
        //        Perform swap provisionally on diagonal block in D.
        //
        Rlarfx("L", 3, 3, u, tau, d, ldd, work);
        Rlarfx("R", 3, 3, u, tau, d, ldd, work);
        //
        //        Test whether to reject swap.
        //
        if (max({abs(d[(2 - 1)]), abs(d[(3 - 1)]), abs(d[(1 - 1)] - t33)}) > thresh) {
            goto statement_50;
        }
        //
        //        Accept swap: apply transformation to the entire matrix T.
        //
        Rlarfx("R", j3, 3, u, tau, &t[(j1 - 1) * ldt], ldt, work);
        Rlarfx("L", 3, n - j1, u, tau, &t[(j1 - 1) + (j2 - 1) * ldt], ldt, work);
        //
        t[(j1 - 1) + (j1 - 1) * ldt] = t33;
        t[(j2 - 1) + (j1 - 1) * ldt] = zero;
        t[(j3 - 1) + (j1 - 1) * ldt] = zero;
        //
        if (wantq) {
            //
            //           Accumulate transformation in the matrix Q.
            //
            Rlarfx("R", n, 3, u, tau, &q[(j1 - 1) * ldq], ldq, work);
        }
        goto statement_40;
    //
    statement_30:
        //
        //        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
        //        that:
        //
        //        H(2) H(1) (  -X11  -X12 ) = (  *  * )
        //                  (  -X21  -X22 )   (  0  * )
        //                  ( scale    0  )   (  0  0 )
        //                  (    0  scale )   (  0  0 )
        //
        u1[1 - 1] = -x[(1 - 1)];
        u1[2 - 1] = -x[(2 - 1)];
        u1[3 - 1] = scale;
        Rlarfg(3, u1[1 - 1], &u1[2 - 1], 1, tau1);
        u1[1 - 1] = one;
        //
        temp = -tau1 * (x[(2 - 1) * ldx] + u1[2 - 1] * x[(2 - 1) + (2 - 1) * ldx]);
        u2[1 - 1] = -temp * u1[2 - 1] - x[(2 - 1) + (2 - 1) * ldx];
        u2[2 - 1] = -temp * u1[3 - 1];
        u2[3 - 1] = scale;
        Rlarfg(3, u2[1 - 1], &u2[2 - 1], 1, tau2);
        u2[1 - 1] = one;
        //
        //        Perform swap provisionally on diagonal block in D.
        //
        Rlarfx("L", 3, 4, u1, tau1, d, ldd, work);
        Rlarfx("R", 4, 3, u1, tau1, d, ldd, work);
        Rlarfx("L", 3, 4, u2, tau2, &d[(2 - 1)], ldd, work);
        Rlarfx("R", 4, 3, u2, tau2, &d[(2 - 1) * ldd], ldd, work);
        //
        //        Test whether to reject swap.
        //
        if (max({abs(d[(3 - 1)]), abs(d[(3 - 1) + (2 - 1) * ldd]), abs(d[(4 - 1)]), abs(d[(4 - 1) + (2 - 1) * ldd])}) > thresh) {
            goto statement_50;
        }
        //
        //        Accept swap: apply transformation to the entire matrix T.
        //
        Rlarfx("L", 3, n - j1 + 1, u1, tau1, &t[(j1 - 1) + (j1 - 1) * ldt], ldt, work);
        Rlarfx("R", j4, 3, u1, tau1, &t[(j1 - 1) * ldt], ldt, work);
        Rlarfx("L", 3, n - j1 + 1, u2, tau2, &t[(j2 - 1) + (j1 - 1) * ldt], ldt, work);
        Rlarfx("R", j4, 3, u2, tau2, &t[(j2 - 1) * ldt], ldt, work);
        //
        t[(j3 - 1) + (j1 - 1) * ldt] = zero;
        t[(j3 - 1) + (j2 - 1) * ldt] = zero;
        t[(j4 - 1) + (j1 - 1) * ldt] = zero;
        t[(j4 - 1) + (j2 - 1) * ldt] = zero;
        //
        if (wantq) {
            //
            //           Accumulate transformation in the matrix Q.
            //
            Rlarfx("R", n, 3, u1, tau1, &q[(j1 - 1) * ldq], ldq, work);
            Rlarfx("R", n, 3, u2, tau2, &q[(j2 - 1) * ldq], ldq, work);
        }
    //
    statement_40:
        //
        if (n2 == 2) {
            //
            //           Standardize new 2-by-2 block T11
            //
            Rlanv2(t[(j1 - 1) + (j1 - 1) * ldt], t[(j1 - 1) + (j2 - 1) * ldt], t[(j2 - 1) + (j1 - 1) * ldt], t[(j2 - 1) + (j2 - 1) * ldt], wr1, wi1, wr2, wi2, cs, sn);
            Rrot(n - j1 - 1, &t[(j1 - 1) + ((j1 + 2) - 1) * ldt], ldt, &t[(j2 - 1) + ((j1 + 2) - 1) * ldt], ldt, cs, sn);
            Rrot(j1 - 1, &t[(j1 - 1) * ldt], 1, &t[(j2 - 1) * ldt], 1, cs, sn);
            if (wantq) {
                Rrot(n, &q[(j1 - 1) * ldq], 1, &q[(j2 - 1) * ldq], 1, cs, sn);
            }
        }
        //
        if (n1 == 2) {
            //
            //           Standardize new 2-by-2 block T22
            //
            j3 = j1 + n2;
            j4 = j3 + 1;
            Rlanv2(t[(j3 - 1) + (j3 - 1) * ldt], t[(j3 - 1) + (j4 - 1) * ldt], t[(j4 - 1) + (j3 - 1) * ldt], t[(j4 - 1) + (j4 - 1) * ldt], wr1, wi1, wr2, wi2, cs, sn);
            if (j3 + 2 <= n) {
                Rrot(n - j3 - 1, &t[(j3 - 1) + ((j3 + 2) - 1) * ldt], ldt, &t[(j4 - 1) + ((j3 + 2) - 1) * ldt], ldt, cs, sn);
            }
            Rrot(j3 - 1, &t[(j3 - 1) * ldt], 1, &t[(j4 - 1) * ldt], 1, cs, sn);
            if (wantq) {
                Rrot(n, &q[(j3 - 1) * ldq], 1, &q[(j4 - 1) * ldq], 1, cs, sn);
            }
        }
        //
    }
    return;
//
//     Exit with INFO = 1 if swap was rejected.
//
statement_50:
    info = 1;
    //
    //     End of Rlaexc
    //
}
