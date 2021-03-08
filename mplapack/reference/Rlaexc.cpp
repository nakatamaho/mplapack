/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaexc.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
#define MTRUE  1

void Rlaexc(INTEGER wantq, INTEGER n, REAL * t, INTEGER ldt, REAL * q, INTEGER ldq, INTEGER j1, INTEGER n1, INTEGER n2, REAL * work, INTEGER * info)
{
    REAL d[16] /* was [4][4] */ ;
    INTEGER k;
    REAL u[3], x[4] /* was [2][2] */ ;
    INTEGER j2, j3, j4;
    REAL u1[3], u2[3];
    INTEGER nd;
    REAL cs, t11, t22, t33, sn, wi1, wi2, wr1, wr2, eps, tau, tau1, tau2;
    INTEGER ierr;
    REAL temp;
    REAL scale, dnorm, xnorm;
    REAL thresh, smlnum;
    REAL One = 1.0, Zero = 0.0, Ten = 10.0;
    REAL mtemp1, mtemp2, mtemp3;

    *info = 0;
//Quick return if possible
    if (n == 0 || n1 == 0 || n2 == 0) {
	return;
    }
    if (j1 + n1 > n) {
	return;
    }
    j2 = j1 + 1;
    j3 = j1 + 2;
    j4 = j1 + 3;
    t11 = 0.0;
    if (n1 == 1 && n2 == 1) {
//Swap two 1-by-1 blocks.
	t11 = t[j1 + j1 * ldt];
	t22 = t[j2 + j2 * ldt];
//Determine the transformation to perform the interchange.
	Rlartg(t[j1 + j2 * ldt], t22 - t11, &cs, &sn, &temp);
//Apply transformation to the matrix T.
	if (j3 <= n) {
	    Rrot(n - j1 - 1, &t[j1 + j3 * ldt], ldt, &t[j2 + j3 * ldt], ldt, cs, sn);
	}
	Rrot(j1 - 1, &t[j1 * ldt + 1], 1, &t[j2 * ldt + 1], 1, cs, sn);
	t[j1 + j1 * ldt] = t22;
	t[j2 + j2 * ldt] = t11;
	if (wantq) {
//Accumulate transformation in the matrix Q.
	    Rrot(n, &q[j1 * ldq + 1], 1, &q[j2 * ldq + 1], 1, cs, sn);
	}
    } else {
//Swapping involves at least one 2-by-2 block.
//Copy the diagonal block of order N1+N2 to the local array D
//and compute its norm.
	nd = n1 + n2;
	Rlacpy("Full", nd, nd, &t[j1 + j1 * ldt], ldt, d, 4);
	dnorm = Rlange("Max", nd, nd, d, 4, &work[0]);
//Compute machine-dependent threshold for test for accepting
//swap.

	eps = Rlamch("P");
	smlnum = Rlamch("S") / eps;
	mtemp1 = eps * Ten * dnorm;
	thresh = max(mtemp1, smlnum);
//Solve T11*X - X*T22 = scale*T12 for X.
	Rlasy2(MFALSE, MFALSE, -1, n1, n2, d, 4, &d[n1 + 1 + (n1 + 1 * 4) - 5], 4, &d[(n1 + 1 * 4) - 4], 4, &scale, x, 2, &xnorm, &ierr);
//Swap the adjacent diagonal blocks.
	k = n1 + n1 + n2 - 3;
	switch (k) {
	case 1:
	    goto L10;
	case 2:
	    goto L20;
	case 3:
	    goto L30;
	}
      L10:
//N1 = 1, N2 = 2: generate elementary reflector H so that:
//( scale, X11, X12 ) H = ( 0, 0, * )
	u[0] = scale;
	u[1] = x[0];
	u[2] = x[2];
	Rlarfg(3, &u[2], u, 1, &tau);
	u[2] = One;
	t11 = t[j1 + j1 * ldt];
//Perform swap provisionally on diagonal block in D.
	Rlarfx("L", 3, 3, u, tau, d, 4, &work[0]);
	Rlarfx("R", 3, 3, u, tau, d, 4, &work[0]);
//Test whether to reject swap.
	mtemp1 = abs(d[2]), mtemp2 = abs(d[6]), mtemp1 = max(mtemp1, mtemp2);
	mtemp2 = abs(d[10] - t11);
	if (max(mtemp1, mtemp2) > thresh) {
	    goto L50;
	}
//Accept swap: apply transformation to the entire matrix T.
	Rlarfx("L", 3, n - j1 + 1, u, tau, &t[j1 + j1 * ldt], ldt, &work[0]);
	Rlarfx("R", j2, 3, u, tau, &t[j1 * ldt + 1], ldt, &work[0]);
	t[j3 + j1 * ldt] = Zero;
	t[j3 + j2 * ldt] = Zero;
	t[j3 + j3 * ldt] = t11;
	if (wantq) {
//Accumulate transformation in the matrix Q.
	    Rlarfx("R", n, 3, u, tau, &q[j1 * ldq + 1], ldq, &work[0]);
	}
	goto L40;
      L20:
//N1 = 2, N2 = 1: generate elementary reflector H so that:
//H (  -X11 ) = ( * )
//  (  -X21 ) = ( 0 )
//  ( scale ) = ( 0 )
	u[0] = -x[0];
	u[1] = -x[0];
	u[2] = scale;
	Rlarfg(3, u, &u[0], 1, &tau);
	u[0] = One;
	t33 = t[j3 + j3 * ldt];
//Perform swap provisionally on diagonal block in D.
	Rlarfx("L", 3, 3, u, tau, d, 4, &work[0]);
	Rlarfx("R", 3, 3, u, tau, d, 4, &work[0]);
//Test whether to reject swap.
	mtemp1 = abs(d[1]), mtemp2 = abs(d[2]), mtemp3 = max(mtemp1, mtemp2);
	mtemp2 = abs(d[0] - t33);
	if (max(mtemp3, mtemp2) > thresh) {
	    goto L50;
	}
//Accept swap: apply transformation to the entire matrix T.
	Rlarfx("R", j3, 3, u, tau, &t[j1 * ldt + 1], ldt, &work[0]);
	Rlarfx("L", 3, n - j1, u, tau, &t[j1 + j2 * ldt], ldt, &work[0]);
	t[j1 + j1 * ldt] = t33;
	t[j2 + j1 * ldt] = Zero;
	t[j3 + j1 * ldt] = Zero;
	if (wantq) {
//Accumulate transformation in the matrix Q.
	    Rlarfx("R", n, 3, u, tau, &q[j1 * ldq + 1], ldq, &work[0]);
	}
	goto L40;

      L30:
//N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
//that:
//H(2) H(1) (  -X11  -X12 ) = (  *  * )
//          (  -X21  -X22 )   (  0  * )
//          ( scale    0  )   (  0  0 )
//          (    0  scale )   (  0  0 )
	u1[0] = -x[0];
	u1[1] = -x[0];
	u1[2] = scale;
	Rlarfg(3, u1, &u1[1], 1, &tau1);
	u1[0] = One;

	temp = -tau1 * (x[2] + u1[1] * x[3]);
	u2[0] = -temp * u1[1] - x[3];
	u2[1] = -temp * u1[2];
	u2[2] = scale;
	Rlarfg(3, u2, &u2[1], 1, &tau2);
	u2[0] = One;
//Perform swap provisionally on diagonal block in D.
	Rlarfx("L", 3, 4, u1, tau1, d, 4, &work[0]);
	Rlarfx("R", 4, 3, u1, tau1, d, 4, &work[0]);
	Rlarfx("L", 3, 4, u2, tau2, &d[0], 4, &work[0]);
	Rlarfx("R", 4, 3, u2, tau2, &d[4], 4, &work[0]);
//Test whether to reject swap.
	mtemp1 = abs(d[2]), mtemp2 = abs(d[6]);
	mtemp3 = max(mtemp1, mtemp2);
	mtemp2 = abs(d[10] - t11);
	if (max(mtemp2, mtemp3) > thresh) {
	    goto L50;
	}
//Accept swap: apply transformation to the entire matrix T.
	Rlarfx("L", 3, n - j1 + 1, u1, tau1, &t[j1 + j1 * ldt], ldt, &work[0]);
	Rlarfx("R", j4, 3, u1, tau1, &t[j1 * ldt + 1], ldt, &work[0]);
	Rlarfx("L", 3, n - j1 + 1, u2, tau2, &t[j2 + j1 * ldt], ldt, &work[0]);
	Rlarfx("R", j4, 3, u2, tau2, &t[j2 * ldt + 1], ldt, &work[0]);

	t[j3 + j1 * ldt] = Zero;
	t[j3 + j2 * ldt] = Zero;
	t[j4 + j1 * ldt] = Zero;
	t[j4 + j2 * ldt] = Zero;
	if (wantq) {
//Accumulate transformation in the matrix Q.
	    Rlarfx("R", n, 3, u1, tau1, &q[j1 * ldq + 1], ldq, &work[0]);
	    Rlarfx("R", n, 3, u2, tau2, &q[j2 * ldq + 1], ldq, &work[0]);
	}
      L40:
	if (n2 == 2) {
//Standardize new 2-by-2 block T11
	    Rlanv2(&t[j1 + j1 * ldt], &t[j1 + j2 * ldt], &t[j2 + j1 * ldt], &t[j2 + j2 * ldt], &wr1, &wi1, &wr2, &wi2, &cs, &sn);
	    Rrot(n - j1 - 1, &t[j1 + (j1 + 2) * ldt], ldt, &t[j2 + (j1 + 2)
							      * ldt], ldt, cs, sn);
	    Rrot(j1 - 1, &t[j1 * ldt + 1], 1, &t[j2 * ldt + 1], 1, cs, sn);
	    if (wantq) {
		Rrot(n, &q[j1 * ldq + 1], 1, &q[j2 * ldq + 1], 1, cs, sn);
	    }
	}
	if (n1 == 2) {
//Standardize new 2-by-2 block T22
	    j3 = j1 + n2;
	    j4 = j3 + 1;
	    Rlanv2(&t[j3 + j3 * ldt], &t[j3 + j4 * ldt], &t[j4 + j3 * ldt], &t[j4 + j4 * ldt], &wr1, &wi1, &wr2, &wi2, &cs, &sn);
	    if (j3 + 2 <= n) {
		Rrot(n - j3 - 1, &t[j3 + (j3 + 2) * ldt], ldt, &t[j4 + (j3 + 2)
								  * ldt], ldt, cs, sn);
	    }
	    Rrot(j3 - 1, &t[j3 * ldt + 1], 1, &t[j4 * ldt + 1], 1, cs, sn);
	    if (wantq) {
		Rrot(n, &q[j3 * ldq + 1], 1, &q[j4 * ldq + 1], 1, cs, sn);
	    }
	}

    }
    return;
//Exit with INFO = 1 if swap was rejected.
  L50:
    *info = 1;
    return;
}
