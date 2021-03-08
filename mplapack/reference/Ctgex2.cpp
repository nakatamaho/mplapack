/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctgex2.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Ctgex2(INTEGER wantq, INTEGER wantz, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, COMPLEX * q, INTEGER ldq, COMPLEX * z, INTEGER ldz, INTEGER j1, INTEGER * info)
{
    COMPLEX f, g;
    INTEGER i, m;
    COMPLEX s[4], t[4];
    REAL cq, sa, sb, cz;
    COMPLEX sq;
    REAL ss, ws;
    COMPLEX sz;
    REAL eps, sum;
    INTEGER weak;
    COMPLEX cdum, work[8];
    REAL scale;
    INTEGER dtrong;
    REAL thresh;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0, Ten = 10.0;
    REAL mtemp1;

    *info = 0;
//Quick return if possible
    if (n <= 1) {
	return;
    }
    m = 2;
    weak = MFALSE;
    dtrong = MFALSE;
//Make a local copy of selected block in (A, B)
    Clacpy("Full", m, m, &A[j1 + j1 * lda], lda, s, 2);
    Clacpy("Full", m, m, &B[j1 + j1 * ldb], ldb, t, 2);
//Compute the threshold for testing the acceptance of swapping.
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    scale = Zero;
    sum = One;
    Clacpy("Full", m, m, s, 2, work, m);
    Clacpy("Full", m, m, t, 2, &work[m * m], m);
    Classq((m * 2) * m, work, 1, &scale, &sum);
    sa = scale * sqrt(sum);

    mtemp1 = eps * Ten * sa;
    thresh = max(mtemp1, smlnum);
//Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
//using Givens rotations and perform the swap tentatively.
    f = s[3] * t[0] - t[3] * s[0];
    g = s[3] * t[2] - t[3] * s[2];
    sa = abs(s[3]);
    sb = abs(t[3]);
    Clartg(g, f, &cz, &sz, &cdum);
    sz = -sz;
    Crot(2, s, 1, &s[2], 1, cz, conj(sz));
    Crot(2, t, 1, &t[2], 1, cz, conj(sz));
    if (sa >= sb) {
	Clartg(s[0], s[1], &cq, &sq, &cdum);
    } else {
	Clartg(t[0], t[1], &cq, &sq, &cdum);
    }
    Crot(2, s, 2, &s[1], 2, cq, sq);
    Crot(2, t, 2, &t[1], 2, cq, sq);
//Weak stability test: |S21| + |T21| <= O(EPS F-norm((S, T)))
    ws = abs(s[1]) + abs(t[1]);
    weak = ws <= thresh;
    if (!weak) {
	goto L20;
    }
    if (MTRUE) {
//Strong stability test:
//   F-norm((A-QL'*S*QR, B-QL'*T*QR)) <= O(EPS*F-norm((A, B)))
	Clacpy("Full", m, m, s, 2, work, m);
	Clacpy("Full", m, m, t, 2, &work[m * m], m);
	Crot(2, work, 1, &work[2], 1, cz, -conj(sz));
	Crot(2, &work[4], 1, &work[6], 1, cz, -conj(sz));
	Crot(2, work, 2, &work[0], 2, cq, -sq);
	Crot(2, &work[4], 2, &work[5], 2, cq, -sq);
	for (i = 0; i < 2; i++) {
	    work[i - 1] = work[i - 1] - A[j1 + i - 1 + j1 * lda];
	    work[i + 1] = work[i + 1] - A[j1 + i - 1 + (j1 + 1) * lda];
	    work[i + 3] = work[i + 3] - B[j1 + i - 1 + j1 * ldb];
	    work[i + 5] = work[i + 5] - B[j1 + i - 1 + (j1 + 1) * ldb];
	}
	scale = Zero;
	sum = One;
	Classq((m * 2) * m, work, 1, &scale, &sum);
	ss = scale * sqrt(sum);
	dtrong = ss <= thresh;
	if (!dtrong) {
	    goto L20;
	}
    }
//If the swap is accepted ("weakly" and "strongly"), apply the
//equivalence transformations to the original matrix pair (A,B)
    Crot(j1 + 1, &A[j1 * lda], 1, &A[(j1 + 1) * lda], 1, cz, conj(sz));
    Crot(j1 + 1, &B[j1 * ldb + 1], 1, &B[(j1 + 1) * ldb + 1], 1, cz, conj(sz));
    Crot(n - j1 + 1, &A[j1 + j1 * lda], lda, &A[j1 + 1 + j1 * lda], lda, cq, sq);
    Crot(n - j1 + 1, &B[j1 + j1 * ldb], ldb, &B[j1 + 1 + j1 * ldb], ldb, cq, sq);
//Set  N1 by N2 (2,1) blocks to 0
    A[j1 + 1 + j1 * lda] = Zero;
    B[j1 + 1 + j1 * ldb] = Zero;
//Accumulate transformations into Q and Z if requested.
    if (wantz) {
	Crot(n, &z[j1 * ldz + 1], 1, &z[(j1 + 1) * ldz + 1], 1, cz, conj(sz));
    }
    if (wantq) {
	Crot(n, &q[j1 * ldq + 1], 1, &q[(j1 + 1) * ldq + 1], 1, cq, conj(sq));
    }
//Exit with INFO = 0 if swap was successfully performed.
    return;
//Exit with INFO = 1 if swap was rejected.
  L20:
    *info = 1;
    return;
}
