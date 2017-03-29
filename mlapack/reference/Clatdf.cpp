/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clatdf.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clatdf(INTEGER ijob, INTEGER n, COMPLEX * z, INTEGER ldz, COMPLEX * rhs, REAL * rdsum, REAL * rdscal, INTEGER * ipiv, INTEGER * jpiv)
{
    INTEGER i, j, k;
    COMPLEX bm, bp, xm[2], xp[2];
    INTEGER info;
    COMPLEX temp, work[8];
    REAL scale;
    COMPLEX pmone;
    REAL rtemp, sminu, rwork[2];
    REAL splus;
    REAL One = 1.0, Zero = 0.0;

    if (ijob != 2) {
	Claswp(1, &rhs[0], ldz, 1, n - 1, &ipiv[0], 1);
//Solve for L-part choosing RHS either to +1 or -One
	pmone = -One;
	for (j = 0; j < n - 1; j++) {
	    bp = rhs[j] + One;
	    bm = rhs[j] - One;
	    splus = One;
//Lockahead for L- part RHS(1:N-1) = +-1
//SPLUS and SMIN computed more efficiently than in BSOLVE[0].
	    splus = splus + Cdotc(n - j, &z[j + 1 + j * ldz], 1, &z[j + 1 + j * ldz], 1).real();
	    sminu = Cdotc(n - j, &z[j + 1 + j * ldz], 1, &rhs[j + 1], 1).real();
	    splus = splus * rhs[j].real();
	    if (splus > sminu) {
		rhs[j] = bp;
	    } else if (sminu > splus) {
		rhs[j] = bm;
	    } else {
//In this case the updating sums are equal and we can
//choose RHS(J) +1 or -One The first time this happens we
//choose -1, thereafter +One This is a simple way to get
//good estimates of matrices like Byers well-known example
//(see [0]). (Not done in BSOLVE.)
		rhs[j] = rhs[j] + pmone;
		pmone = One;
	    }
//Compute the remaining r.h.s.
	    temp = -rhs[j];
	    Caxpy(n - j, temp, &z[j + 1 + j * ldz], 1, &rhs[j + 1], 1);
	}
//Solve for U- part, lockahead for RHS(N) = +-One This is not done
//In BSOLVE and will hopefully give us a better estimate because
//any ill-conditioning of the original matrix is transfered to U
//and not to L. U(N, N) is an approximation to sigma_min(LU).
	Ccopy(n - 1, &rhs[0], 1, work, 1);
	work[n] = rhs[n] + One;
	rhs[n] = rhs[n] - One;
	splus = Zero;
	sminu = Zero;
	for (i = n; i >= 1; i--) {
	    temp = One / z[i + i * ldz];
	    work[i] = work[i] * temp;
	    rhs[i] = rhs[i] * temp;
	    for (k = i + 1; k <= n; k++) {
		work[i] = work[i] - work[k] * (z[i + k * ldz] * temp);
		rhs[i] = rhs[i] - rhs[k] * (z[i + k * ldz] * temp);
	    }
	    splus = splus + abs(work[i - 1]);
	    sminu = sminu + abs(rhs[i]);
	}
	if (splus > sminu) {
	    Ccopy(n, work, 1, &rhs[0], 1);
	}
//Apply the permutations JPIV to the computed solution (RHS)
	Claswp(1, &rhs[0], ldz, 1, n - 1, &jpiv[0], -1);
//Compute the sum of squares
	Classq(n, &rhs[0], 1, rdscal, rdsum);
	return;
    }
//ENTRY IJOB = 2
//Compute approximate nullvector XM of Z
    Cgecon("I", n, &z[0], ldz, One, &rtemp, &work[0], rwork, &info);
    Ccopy(n, &work[n], 1, xm, 1);
//Compute RHS
    Claswp(1, xm, ldz, 1, n - 1, &ipiv[0], -1);
    temp = One / sqrt(Cdotc(n, xm, 1, xm, 1));
    Cscal(n, temp, xm, 1);
    Ccopy(n, xm, 1, xp, 1);
    Caxpy(n, (COMPLEX) One, &rhs[0], 1, xp, 1);
    Caxpy(n, -(COMPLEX) One, xm, 1, &rhs[0], 1);
    Cgesc2(n, &z[0], ldz, &rhs[0], &ipiv[0], &jpiv[0], &scale);
    Cgesc2(n, &z[0], ldz, xp, &ipiv[0], &jpiv[0], &scale);
    if (RCasum(n, xp, 1) > RCasum(n, &rhs[0], 1)) {
	Ccopy(n, xp, 1, &rhs[0], 1);
    }
//Compute the sum of squares
    Classq(n, &rhs[0], 1, rdscal, rdsum);
    return;
}
