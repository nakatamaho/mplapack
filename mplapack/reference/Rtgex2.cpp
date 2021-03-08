/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgex2.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mpblas.h>
#include <mplapack.h>

#define MFALSE 0
#define MTRUE 1

void
Rtgex2(INTEGER wantq, INTEGER wantz, INTEGER n, REAL * A, INTEGER lda, REAL * B,
       INTEGER ldb, REAL * q, INTEGER ldq, REAL * z, INTEGER ldz, INTEGER j1, INTEGER n1, INTEGER n2, REAL * work, INTEGER lwork, INTEGER * info)
{
    REAL f, g;
    INTEGER i, m;
    REAL s[16], t[16], be[2], ai[2], ar[2], sa, sb, li[16], ir[16], ss, ws, eps;
    INTEGER weak;
    REAL ddum;
    INTEGER idum;
    REAL taul[4], dsum;
    REAL taur[4], scpy[16], tcpy[16];
    REAL scale, bqra21, brqa21;
    REAL licop[16];
    INTEGER linfo;
    REAL ircop[16], dnorm;
    INTEGER iwork[4];
    REAL dscale;
    INTEGER dtrong;
    REAL thresh, smlnum;
    REAL One = 1.0, Zero = 0.0, Ten = 10.0;
    REAL mtemp1, mtemp2;

    *info = 0;
//Quick return if possible
    if (n <= 1 || n1 <= 0 || n2 <= 0) {
	return;
    }
    if (n1 > n || j1 + n1 > n) {
	return;
    }
    m = n1 + n2;
    if (lwork < max(max((INTEGER) 1, n * m), m * m << 1)) {
	*info = -16;
	work[1] = max(max((INTEGER) 1, n * m), m * m << 1);
	return;
    }

    weak = MFALSE;
    dtrong = MFALSE;
//Make a local copy of selected block
    Rlaset("Full", 4, 4, Zero, Zero, li, 4);
    Rlaset("Full", 4, 4, Zero, Zero, ir, 4);
    Rlacpy("Full", m, m, &A[j1 + j1 * lda], lda, s, 4);
    Rlacpy("Full", m, m, &B[j1 + j1 * ldb], ldb, t, 4);
//Compute threshold for testint acceptance of swappint.
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    dscale = Zero;
    dsum = One;
    Rlacpy("Full", m, m, s, 4, &work[0], m);
    Rlassq(m * m, &work[0], 1, &dscale, &dsum);
    Rlacpy("Full", m, m, t, 4, &work[0], m);
    Rlassq(m * m, &work[0], 1, &dscale, &dsum);
    dnorm = dscale * sqrt(dsum);
    mtemp1 = eps * Ten * dnorm;
    mtemp2 = smlnum;
    thresh = max(mtemp1, mtemp2);
    if (m == 2) {
//CASE 1: Swap 1-by-1 and 1-by-1 blocks.
//Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
//usint Givens rotations and perform the swap tentatively.
	f = s[5] * t[0] - t[5] * s[0];
	g = s[5] * t[4] - t[5] * s[4];
	sb = abs(t[5]);
	sa = abs(s[5]);
	Rlartg(f, g, &ir[4], ir, &ddum);
	ir[1] = -ir[4];
	ir[5] = ir[0];
	Rrot(2, s, 1, &s[4], 1, ir[0], ir[1]);
	Rrot(2, t, 1, &t[4], 1, ir[0], ir[1]);
	if (sa >= sb) {
	    Rlartg(s[0], s[1], &li[0], &li[1], &ddum);
	} else {
	    Rlartg(t[0], t[1], &li[0], &li[1], &ddum);
	}
	Rrot(2, s, 4, &s[1], 4, li[0], li[1]);
	Rrot(2, t, 4, &t[1], 4, li[0], li[1]);
	li[5] = li[0];
	li[4] = -li[1];
//Weak stability test:
//|S21| + |T21| <= O(EPS * F-norm((S, T)))
	ws = abs(s[1]) + abs(t[1]);
	weak = ws <= thresh;
	if (!weak) {
	    goto L70;
	}
	if (MTRUE) {
//Strong stability test:
//  F-norm((A-QL'*S*QR, B-QL'*T*QR)) <= O(EPS*F-norm((A,B)))
	    Rlacpy("Full", m, m, &A[j1 + j1 * lda], lda, &work[m * m + 1], m);
	    Rgemm("N", "N", m, m, m, One, li, 4, s, 4, Zero, &work[0], m);
	    Rgemm("N", "T", m, m, m, -One, &work[0], m, ir, 4, One, &work[m * m + 1], m);
	    dscale = Zero;
	    dsum = One;
	    Rlassq(m * m, &work[m * m + 1], 1, &dscale, &dsum);
	    Rlacpy("Full", m, m, &B[j1 + j1 * ldb], ldb, &work[m * m + 1], m);
	    Rgemm("N", "N", m, m, m, One, li, 4, t, 4, Zero, &work[0], m);
	    Rgemm("N", "T", m, m, m, -One, &work[0], m, ir, 4, One, &work[m * m + 1], m);
	    Rlassq(m * m, &work[m * m + 1], 1, &dscale, &dsum);
	    ss = dscale * sqrt(dsum);
	    dtrong = ss <= thresh;
	    if (!dtrong) {
		goto L70;
	    }
	}
//Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
//       (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
	Rrot(j1 + 1, &A[j1 * lda], 1, &A[(j1 + 1) * lda], 1, ir[0], ir[1]);
	Rrot(j1 + 1, &B[j1 * ldb + 1], 1, &B[(j1 + 1) * ldb + 1], 1, ir[0], ir[1]);
	Rrot(n - j1 + 1, &A[j1 + j1 * lda], lda, &A[j1 + 1 + j1 * lda], lda, li[0], li[1]);
	Rrot(n - j1 + 1, &B[j1 + j1 * ldb], ldb, &B[j1 + 1 + j1 * ldb], ldb, li[0], li[1]);
//Set  N1-by-N2 (2,1) - blocks to ZERO.
	A[j1 + 1 + j1 * lda] = Zero;
	B[j1 + 1 + j1 * ldb] = Zero;
//Accumulate transformations into Q and Z if requested.
	if (wantz) {
	    Rrot(n, &z[j1 * ldz + 1], 1, &z[(j1 + 1) * ldz + 1], 1, ir[0], ir[1]);
	}
	if (wantq) {
	    Rrot(n, &q[j1 * ldq + 1], 1, &q[(j1 + 1) * ldq + 1], 1, li[0], li[1]);
	}
//Exit with INFO = 0 if swap was successfully performed.
	return;
    } else {
//CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
//        and 2-by-2 blocks.
//Solve the generalized Sylvester equation
//         S11 * R - L * S22 = SCALE * S12
//         T11 * R - L * T22 = SCALE * T12
//for R and L. Solutions in LI and IR.
	Rlacpy("Full", n1, n2, &t[(n1 + 1 * 4) - 4], 4, li, 4);
	Rlacpy("Full", n1, n2, &s[(n1 + 1 * 4) - 4], 4, &ir[n2 + 1 + (n1 + 1 * 4) - 5], 4);
	Rtgsy2("N", 0, n1, n2, s, 4, &s[n1 + 1 + (n1 + 1 * 4) - 5]
	       , 4, &ir[n2 + 1 + (n1 + 1 * 4) - 5], 4, t, 4, &t[n1 + 1 + (n1 + 1 * 4) - 5], 4, li, 4, &scale, &dsum, &dscale, iwork, &idum, &linfo);
//Compute orthogonal matrix QL:
//            QL' * LI = [ TL ]
//                       [ 0  ]
//where
//            LI =  [      -L              ]
//                  [ SCALE * identity(N2) ]
	for (i = 0; i < n2; i++) {
	    Rscal(n1, -One, &li[(i * 4) - 4], 1);
	    li[n1 + i + (i * 4) - 5] = scale;
	}
	Rgeqr2(m, n2, li, 4, taul, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
	Rorg2r(m, m, n2, li, 4, taul, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
//Compute orthogonal matrix RQ:
//            IR * RQ' =   [ 0  TR],
// where IR = [ SCALE * identity(N1), R ]
	for (i = 0; i < n1; i++) {
	    ir[n2 + i + (i * 4) - 5] = scale;
	}
	Rgerq2(n1, m, &ir[n2], 4, taur, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
	Rorgr2(m, m, n1, ir, 4, taur, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
//Perform the swappint tentatively:
	Rgemm("T", "N", m, m, m, One, li, 4, s, 4, Zero, &work[0], m);
	Rgemm("N", "T", m, m, m, One, &work[0], m, ir, 4, Zero, s, 4);
	Rgemm("T", "N", m, m, m, One, li, 4, t, 4, Zero, &work[0], m);
	Rgemm("N", "T", m, m, m, One, &work[0], m, ir, 4, Zero, t, 4);
	Rlacpy("F", m, m, s, 4, scpy, 4);
	Rlacpy("F", m, m, t, 4, tcpy, 4);
	Rlacpy("F", m, m, ir, 4, ircop, 4);
	Rlacpy("F", m, m, li, 4, licop, 4);
//Triangularize the B-part by an RQ factorization.
//Apply transformation (from left) to A-part, givint S.
	Rgerq2(m, m, t, 4, taur, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
	Rormr2("R", "T", m, m, m, t, 4, taur, s, 4, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
	Rormr2("L", "N", m, m, m, t, 4, taur, ir, 4, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
//Compute F-norm(S21) in BRQA2One (T21 is Zero)
	dscale = Zero;
	dsum = One;
	for (i = 0; i < n2; i++) {
	    Rlassq(n1, &s[n2 + 1 + (i * 4) - 5], 1, &dscale, &dsum);
	}
	brqa21 = dscale * sqrt(dsum);
//Triangularize the B-part by a QR factorization.
//Apply transformation (from right) to A-part, givint S.
	Rgeqr2(m, m, tcpy, 4, taul, &work[0], &linfo);
	if (linfo != 0) {
	    goto L70;
	}
	Rorm2r("L", "T", m, m, m, tcpy, 4, taul, scpy, 4, &work[0], info);
	Rorm2r("R", "N", m, m, m, tcpy, 4, taul, licop, 4, &work[0], info);
	if (linfo != 0) {
	    goto L70;
	}
//Compute F-norm(S21) in BQRA2One (T21 is Zero)
	dscale = Zero;
	dsum = One;
	for (i = 0; i < n2; i++) {
	    Rlassq(n1, &scpy[n2 + 1 + (i * 4) - 5], 1, &dscale, &dsum);
	}
	bqra21 = dscale * sqrt(dsum);
//Decide which method to use.
//  Weak stability test:
//     F-norm(S21) <= O(EPS * F-norm((S, T)))
	if (bqra21 <= brqa21 && bqra21 <= thresh) {
	    Rlacpy("F", m, m, scpy, 4, s, 4);
	    Rlacpy("F", m, m, tcpy, 4, t, 4);
	    Rlacpy("F", m, m, ircop, 4, ir, 4);
	    Rlacpy("F", m, m, licop, 4, li, 4);
	} else if (brqa21 >= thresh) {
	    goto L70;
	}
//Set lower triangle of B-part to zero
	Rlaset("Lower", m - 1, m - 1, Zero, Zero, &t[1], 4);
	if (MTRUE) {
//Strong stability test:
//   F-norm((A-QL*S*QR', B-QL*T*QR')) <= O(EPS*F-norm((A,B)))
	    Rlacpy("Full", m, m, &A[j1 + j1 * lda], lda, &work[m * m + 1], m);
	    Rgemm("N", "N", m, m, m, One, li, 4, s, 4, Zero, &work[0], m);
	    Rgemm("N", "N", m, m, m, -One, &work[0], m, ir, 4, One, &work[m * m + 1], m);
	    dscale = Zero;
	    dsum = One;
	    Rlassq(m * m, &work[m * m + 1], 1, &dscale, &dsum);
	    Rlacpy("Full", m, m, &B[j1 + j1 * ldb], ldb, &work[m * m + 1], m);
	    Rgemm("N", "N", m, m, m, One, li, 4, t, 4, Zero, &work[0], m);
	    Rgemm("N", "N", m, m, m, -One, &work[0], m, ir, 4, One, &work[m * m + 1], m);
	    Rlassq(m * m, &work[m * m + 1], 1, &dscale, &dsum);
	    ss = dscale * sqrt(dsum);
	    dtrong = ss <= thresh;
	    if (!dtrong) {
		goto L70;
	    }
	}
//If the swap is accepted ("weakly" and "strongly"), apply the
//transformations and set N1-by-N2 (2,1)-block to zero.
	Rlaset("Full", n1, n2, Zero, Zero, &s[n2], 4);
//copy back M-by-M diagonal block startint at index J1 of (A, B)
	Rlacpy("F", m, m, s, 4, &A[j1 + j1 * lda], lda);
	Rlacpy("F", m, m, t, 4, &B[j1 + j1 * ldb], ldb);
	Rlaset("Full", 4, 4, Zero, Zero, t, 4);
//Standardize existing 2-by-2 blocks.
	for (i = 0; i < m * m; i++) {
	    work[i] = Zero;
	}
	work[1] = One;
	t[0] = One;
	idum = lwork - m * m - 2;
	if (n2 > 1) {
	    Rlagv2(&A[j1 + j1 * lda], lda, &B[j1 + j1 * ldb], ldb, ar, ai, be, &work[0], &work[2], t, &t[1]);
	    work[m + 1] = -work[2];
	    work[m + 2] = work[1];
	    t[n2 + (n2 * 4) - 5] = t[0];
	    t[4] = -t[1];
	}
	work[m * m] = One;
	t[m + (m * 4) - 5] = One;

	if (n1 > 1) {
	    Rlagv2(&A[j1 + n2 + (j1 + n2) * lda], lda, &B[j1 + n2 +
							  (j1 + n2) * ldb], ldb, taur, taul, &work[m * m + 1],
		   &work[n2 * m + n2 + 1], &work[n2 * m + n2 + 2], &t[n2 + 1 + (n2 + 1 * 4) - 5], &t[m + (m - 1 * 4) - 5]);
	    work[m * m] = work[n2 * m + n2 + 1];
	    work[m * m - 1] = -work[n2 * m + n2 + 2];
	    t[m + (m * 4) - 5] = t[n2 + 1 + (n2 + 1 * 4) - 5];
	    t[m - 1 + (m * 4) - 5] = -t[m + (m - 1 * 4) - 5];
	}
	Rgemm("T", "N", n2, n1, n2, One, &work[0], m, &A[j1 + (j1 + n2) * lda], lda, Zero, &work[m * m + 1], n2);
	Rlacpy("Full", n2, n1, &work[m * m + 1], n2, &A[j1 + (j1 + n2) * lda], lda);
	Rgemm("T", "N", n2, n1, n2, One, &work[0], m, &B[j1 + (j1 + n2) * ldb], ldb, Zero, &work[m * m + 1], n2);
	Rlacpy("Full", n2, n1, &work[m * m + 1], n2, &B[j1 + (j1 + n2) * ldb], ldb);
	Rgemm("N", "N", m, m, m, One, li, 4, &work[0], m, Zero, &work[m * m + 1], m);
	Rlacpy("Full", m, m, &work[m * m + 1], m, li, 4);
	Rgemm("N", "N", n2, n1, n1, One, &A[j1 + (j1 + n2) * lda], lda, &t[n2 + 1 + (n2 + 1 * 4) - 5], 4, Zero, &work[0], n2);
	Rlacpy("Full", n2, n1, &work[0], n2, &A[j1 + (j1 + n2) * lda], lda);
	Rgemm("N", "N", n2, n1, n1, One, &B[j1 + (j1 + n2) * ldb], ldb, &t[n2 + 1 + (n2 + 1 * 4) - 5], 4, Zero, &work[0], n2);
	Rlacpy("Full", n2, n1, &work[0], n2, &B[j1 + (j1 + n2) * ldb], ldb);
	Rgemm("T", "N", m, m, m, One, ir, 4, t, 4, Zero, &work[0], m);
	Rlacpy("Full", m, m, &work[0], m, ir, 4);
/*        Accumulate transformations into Q and Z if requested. */
	if (wantq) {
	    Rgemm("N", "N", n, m, m, One, &q[j1 * ldq + 1], ldq, li, 4, Zero, &work[0], n);
	    Rlacpy("Full", n, m, &work[0], n, &q[j1 * ldq + 1], ldq);

	}
	if (wantz) {
	    Rgemm("N", "N", n, m, m, One, &z[j1 * ldz + 1], ldz, ir, 4, Zero, &work[0], n);
	    Rlacpy("Full", n, m, &work[0], n, &z[j1 * ldz + 1], ldz);

	}
//Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and
//       (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)).
	i = j1 + m;
	if (i <= n) {
	    Rgemm("T", "N", m, n - i + 1, m, One, li, 4, &A[j1 + i * lda], lda, Zero, &work[0], m);
	    Rlacpy("Full", m, n - i + 1, &work[0], m, &A[j1 + i * lda], lda);
	    Rgemm("T", "N", m, n - i + 1, m, One, li, 4, &B[j1 + i * ldb], lda, Zero, &work[0], m);
	    Rlacpy("Full", m, n - i + 1, &work[0], m, &B[j1 + i * ldb], ldb);
	}
	if (i > 0) {
	    Rgemm("N", "N", j1 - 1, m, m, One, &A[j1 * lda], lda, ir, 4, Zero, &work[0], i);
	    Rlacpy("Full", j1 - 1, m, &work[0], i, &A[j1 * lda], lda);
	    Rgemm("N", "N", j1 - 1, m, m, One, &B[j1 * ldb + 1], ldb, ir, 4, Zero, &work[0], i);
	    Rlacpy("Full", j1 - 1, m, &work[0], i, &B[j1 * ldb + 1], ldb);
	}
//Exit with INFO = 0 if swap was successfully performed.
	return;
    }
//Exit with INFO = 1 if swap was rejected.
  L70:
    *info = 1;
    return;
}
