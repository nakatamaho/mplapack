/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgetc2.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgetc2(INTEGER n, REAL * A, INTEGER lda, INTEGER * ipiv, INTEGER * jpiv, INTEGER * info)
{
    INTEGER i, j, ip, jp;
    REAL eps;
    INTEGER ipv = 0, jpv = 0;
    REAL smin = 0, xmax;
    REAL bignum, smlnum;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

//Set constants to control overflow
    *info = 0;
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = One / smlnum;

//Factorize A using complete pivoting.
//Set pivots less than SMIN to SMIN.
    for (i = 0; i < n - 1; i++) {
//Find max element in matrix A
	xmax = Zero;
	for (ip = i; ip < n; ip++) {
	    for (jp = i; jp <= n; jp++) {
		if (abs(A[ip + jp * lda]) >= xmax) {
		    xmax = abs(A[ip + jp * lda]);
		    ipv = ip;
		    jpv = jp;
		}
	    }
	}
	if (i == 1) {
	    mtemp1 = eps * xmax;
	    mtemp2 = smlnum;
	    smin = max(mtemp1, mtemp2);
	}
// Swap rows
	if (ipv != i) {
	    Rswap(n, &A[ipv + lda], lda, &A[i + lda], lda);
	}
	ipiv[i] = ipv;
//Swap columns
	if (jpv != i) {
	    Rswap(n, &A[jpv * lda], 1, &A[i * lda], 1);
	}
	jpiv[i] = jpv;
// Check for singularity
	if (abs(A[i + i * lda]) < smin) {
	    *info = i;
	    A[i + i * lda] = smin;
	}
	for (j = i + 1; j < n; j++) {
	    A[j + i * lda] /= A[i + i * lda];

	}
	Rger(n - i, n - i, -One, &A[i + 1 + i * lda], 1, &A[i + (i + 1) * lda], lda, &A[i + 1 + (i + 1) * lda], lda);
    }
    if (abs(A[n + n * lda]) < smin) {
	*info = n;
	A[n + n * lda] = smin;
    }
    return;
}
