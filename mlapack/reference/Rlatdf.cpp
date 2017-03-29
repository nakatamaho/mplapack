/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlatdf.cpp,v 1.6 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>

void Rlatdf(INTEGER ijob, INTEGER n, REAL * z, INTEGER ldz, REAL * rhs, REAL * rdsum, REAL * rdscal, INTEGER * ipiv, INTEGER * jpiv)
{
    INTEGER i, j, k;
    REAL bm, bp, xm[8], xp[8];
    INTEGER info;
    REAL temp, work[32];
    REAL pmone;
    REAL sminu;
    INTEGER iwork[8];
    REAL splus;
    REAL One = 1.0, Zero = 0.0;

    if (ijob != 2) {
//Apply permutations IPIV to RHS
	Rlaswp(1, &rhs[1], ldz, 1, n - 1, &ipiv[1], 1);
//Solve for L-part choosing RHS either to +1 or -One
	pmone = -One;
	for (j = 0; j < n - 1; j++) {
	    bp = rhs[j] + One;
	    bm = rhs[j] - One;
	    splus = One;
//Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and
//SMIN computed more efficiently than in BSOLVE [1].
	    splus = splus + Rdot(n - j, &z[j + 1 + j * ldz], 1, &z[j + 1 + j * ldz], 1);
	    sminu = Rdot(n - j, &z[j + 1 + j * ldz], 1, &rhs[j + 1], 1);
	    splus = splus * rhs[j];
	    if (splus > sminu) {
		rhs[j] = bp;
	    } else if (sminu > splus) {
		rhs[j] = bm;
	    } else {
//In this case the updating sums are equal and we can
//choose RHS(J) +1 or -One The first time this happens
//we choose -1, thereafter +One This is a simple way to
//get good estimates of matrices like Byers well-known
//example (see [1]). (Not done in BSOLVE.)
		rhs[j] = rhs[j] + pmone;
		pmone = One;
	    }
//Compute the remaining r.h.s.
	    Raxpy(n - j, -rhs[j], &z[j + 1 + j * ldz], 1, &rhs[j + 1], 1);
	}
//Solve for U-part, look-ahead for RHS(N) = +-One This is not done
//in BSOLVE and will hopefully give us a better estimate because
//any ill-conditioning of the original matrix is transfered to U
//and not to L. U(N, N) is an approximation to sigma_min(LU).
	Rcopy(n - 1, &rhs[1], 1, xp, 1);
	xp[n - 1] = rhs[n] + One;
	rhs[n] = rhs[n] - One;
	splus = Zero;
	sminu = Zero;
	for (i = n; i >= 1; i--) {
	    temp = One / z[i + i * ldz];
	    xp[i - 1] = xp[i - 1] * temp;
	    rhs[i] = rhs[i] * temp;
	    for (k = i + 1; k <= n; k++) {
		xp[i - 1] = xp[i - 1] - xp[k - 1] * (z[i + k * ldz] * temp);
		rhs[i] = rhs[i] - rhs[k] * (z[i + k * ldz] * temp);
	    }
	    splus = splus + abs(xp[i - 1]);
	    sminu = sminu + abs(rhs[i]);
	}
	if (splus > sminu) {
	    Rcopy(n, xp, 1, &rhs[1], 1);
	}
//Apply the permutations JPIV to the computed solution (RHS)
	Rlaswp(1, &rhs[1], ldz, 1, n - 1, &jpiv[1], -1);
//Compute the sum of squares
	Rlassq(n, &rhs[1], 1, rdscal, rdsum);
    } else {
//IJOB = 2, Compute approximate nullvector XM of Z
	Rgecon("I", n, &z[0], ldz, One, &temp, work, iwork, &info);
	Rcopy(n, &work[n], 1, xm, 1);
//Compute RHS
	Rlaswp(1, xm, ldz, 1, n - 1, &ipiv[1], -1);
	temp = One / sqrt(Rdot(n, xm, 1, xm, 1));
	Rscal(n, temp, xm, 1);
	Rcopy(n, xm, 1, xp, 1);
	Raxpy(n, One, &rhs[1], 1, xp, 1);
	Raxpy(n, -One, xm, 1, &rhs[1], 1);
	Rgesc2(n, &z[0], ldz, &rhs[1], &ipiv[1], &jpiv[1], &temp);
	Rgesc2(n, &z[0], ldz, xp, &ipiv[1], &jpiv[1], &temp);
	if (Rasum(n, xp, 1) > Rasum(n, &rhs[1], 1)) {
	    Rcopy(n, xp, 1, &rhs[1], 1);
	}
//Compute the sum of squares
	Rlassq(n, &rhs[1], 1, rdscal, rdsum);
    }
    return;
}
