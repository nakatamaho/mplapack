/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgeqpf.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgeqpf(INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * jpvt, COMPLEX * tau, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j, ma, mn;
    COMPLEX aii;
    INTEGER pvt;
    REAL temp, temp2, tol3z;
    INTEGER itemp;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;
//Test the input arguments
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Cgeqpf", -(*info));
	return;
    }
    mn = min(m, n);
    tol3z = sqrt(Rlamch("Epsilon"));
//Move initial columns up front
    itemp = 1;
    for (i = 0; i < n; i++) {
	if (jpvt[i] != 0) {
	    if (i != itemp) {
		Cswap(m, &A[i * lda], 1, &A[itemp * lda], 1);
		jpvt[i] = jpvt[itemp];
		jpvt[itemp] = i;
	    } else {
		jpvt[i] = i;
	    }
	    itemp++;
	} else {
	    jpvt[i] = i;
	}
    }
    itemp--;
//Compute the QR factorization and update remaining columns
    if (itemp > 0) {
	ma = min(itemp, m);
	Cgeqr2(m, ma, &A[0], lda, &tau[1], &work[0], info);
	if (ma < n) {
	    Cunm2r("Left", "Conjugate transpose", m, n - ma, ma, &A[0]
		   , lda, &tau[1], &A[(ma + 1) * lda], lda, &work[0], info);
	}
    }
    if (itemp < mn) {
//Initialize partial column norms. The first n elements of
//work store the exact column norms.
	for (i = itemp + 1; i <= n; i++) {
	    rwork[i] = RCnrm2(m - itemp, &A[itemp + 1 + i * lda], 1);
	    rwork[n + i] = rwork[i];
	}
//Compute factorization
	for (i = itemp + 1; i <= mn; i++) {
//Determine ith pivot column and swap if necessary
	    pvt = i - 1 + iRamax(n - i + 1, &rwork[i], 1);
	    if (pvt != i) {
		Cswap(m, &A[pvt * lda], 1, &A[i * lda], 1);
		itemp = jpvt[pvt];
		jpvt[pvt] = jpvt[i];
		jpvt[i] = itemp;
		rwork[pvt] = rwork[i];
		rwork[n + pvt] = rwork[n + i];
	    }
//Generate elementary reflector H(i)
	    aii = A[i + i * lda];
	    Clarfg(m - i + 1, &aii, &A[min(i + 1, m) + i * lda], 1, &tau[i]);
	    A[i + i * lda] = aii;
	    if (i < n) {
//Apply H(i) to A(i:m,i+1:n) from the left
		aii = A[i + i * lda];
		A[i + i * lda] = One;
		Clarf("Left", m - i + 1, n - i, &A[i + i * lda], 1, conj(tau[i]), &A[i + (i + 1) * lda], lda, &work[0]);
		A[i + i * lda] = aii;
	    }
//Update partial column norms
	    for (j = i + 1; j <= n; j++) {
		if (rwork[j] != Zero) {
//NOTE: The following 4 lines follow from the analysis in
//Lapack Working Note 176.
		    temp = abs(A[i + j * lda]) / rwork[j];
		    mtemp1 = Zero, mtemp2 = (temp + One) * (One - temp);
		    temp = max(mtemp1, mtemp2);
		    temp2 = temp * (rwork[j] / rwork[n + j]) * (rwork[j] / rwork[n + j]);
		    if (temp2 <= tol3z) {
			if (m - i > 0) {
			    rwork[j] = RCnrm2(m - i, &A[i + 1 + j * lda], 1);
			    rwork[n + j] = rwork[j];
			} else {
			    rwork[j] = Zero;
			    rwork[n + j] = Zero;
			}
		    } else {
			rwork[j] *= sqrt(temp);
		    }
		}
	    }
	}
    }
    return;
}
