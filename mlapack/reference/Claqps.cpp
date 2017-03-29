/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claqps.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Claqps(INTEGER m, INTEGER n, INTEGER offset, INTEGER nb, INTEGER * kb, COMPLEX * A, INTEGER lda,
       INTEGER * jpvt, COMPLEX * tau, REAL * vn1, REAL * vn2, COMPLEX * auxv, COMPLEX * f, INTEGER ldf)
{
    INTEGER j, k, rk;
    COMPLEX akk;
    INTEGER pvt;
    REAL temp;
    REAL temp2, tol3z;
    INTEGER itemp;
    INTEGER lsticc, lastrk;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

    lastrk = min(m, n + offset);
    lsticc = 0;
    k = 0;
    tol3z = sqrt(Rlamch("Epsilon"));

// Beginning of while loop.
  L10:
    if (k < nb && lsticc == 0) {
	k++;
	rk = offset + k;
//Determine ith pivot column and swap if necessary
	pvt = k - 1 + iRamax(n - k + 1, &vn1[k], 1);
	if (pvt != k) {
	    Cswap(m, &A[pvt * lda], 1, &A[k * lda], 1);
	    Cswap(k - 1, &f[pvt + ldf], ldf, &f[k + ldf], ldf);
	    itemp = jpvt[pvt];
	    jpvt[pvt] = jpvt[k];
	    jpvt[k] = itemp;
	    vn1[pvt] = vn1[k];
	    vn2[pvt] = vn2[k];
	}
//Apply previous Householder reflectors to column K: 
//A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
	if (k > 1) {
	    Cgemv("No transpose", m - rk + 1, k - 1, (COMPLEX) - One, &A[rk + lda], lda, &f[k + ldf], ldf, (COMPLEX) One, &A[rk + k * lda], 1);
	}
//Generate elementary reflector H(k).
	if (rk < m) {
	    Clarfg(m - rk + 1, &A[rk + k * lda], &A[rk + 1 + k * lda], 1, &tau[k]);
	} else {
	    Clarfg(1, &A[rk + k * lda], &A[rk + k * lda], 1, &tau[k]);
	}

	akk = A[rk + k * lda];
	A[rk + k * lda] = One;
//Compute Kth column of F:
//Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K).
	if (k < n) {
	    Cgemv("Transpose", m - rk + 1, n - k, tau[k], &A[rk + (k + 1) * lda], lda, &A[rk + k * lda], 1, Zero, &f[k + 1 + k * ldf], 1);
	}
//Padding F(1:K,K) with zeros.
	for (j = 0; j < k; j++) {
	    f[j + k * ldf] = Zero;
	}
//Incremental updating of F:
//F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'
//            *A(RK:M,K).
	if (k > 1) {
	    Cgemv("Transpose", m - rk + 1, k - 1, -tau[k], &A[rk + lda], lda, &A[rk + k * lda], 1, Zero, &auxv[1], 1);
	    Cgemv("No transpose", n, k - 1, One, &f[ldf + 1], ldf, &auxv[1], 1, One, &f[k * ldf + 1], 1);
	}
//Update the current row of A:
//A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.
	if (k < n) {
	    Cgemv("No transpose", n - k, k, (COMPLEX) - One, &f[k + 1 + ldf], ldf, &A[rk + lda], lda, (COMPLEX) One, &A[rk + (k + 1) * lda], lda);
	}
//Update partial column norms.
	if (rk < lastrk) {
	    for (j = k + 1; j <= n; j++) {
		if (vn1[j] != Zero) {
//NOTE: The following 4 lines follow from the analysis in
//Lapack Working Note 176.
		    temp = abs(A[rk + j * lda]) / vn1[j];
		    mtemp1 = Zero, mtemp2 = (temp + One) * (One - temp);
		    temp = max(mtemp1, mtemp2);
		    temp2 = temp * ((vn1[j] / vn2[j]) * (vn1[j] / vn2[j]));
		    if (temp2 <= tol3z) {
			vn2[j] = lsticc;
			lsticc = j;
		    } else {
			vn1[j] = vn1[j] * sqrt(temp);
		    }
		}

	    }
	}

	A[rk + k * lda] = akk;
//End of while loop.
	goto L10;
    }
    *kb = k;
    rk = offset + *kb;
//Apply the block reflector to the rest of the matrix:
//A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
//                    A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'.
    if (*kb < min(n, m - offset)) {
	Cgemm("No transpose", "Transpose", m - rk, n - *kb, *kb, (COMPLEX) - One, &A[rk + 1 + lda], lda, &f[*kb + 1 + ldf], ldf, One, &A[rk + 1 + (*kb + 1) * lda], lda);
    }
//Recomputation of difficult columns.
  L40:
    if (lsticc > 0) {
	itemp = nint(vn2[lsticc]);
	vn1[lsticc] = RCnrm2(m - rk, &A[rk + 1 + lsticc * lda], 1);
//NOTE: The computation of VN1( LSTICC ) relies on the fact that
//SNRM2 does not fail on vectors with norm below the value of
//SQRT(DLAMCH('S'))
	vn2[lsticc] = vn1[lsticc];
	lsticc = itemp;
	goto L40;
    }
    return;
}
