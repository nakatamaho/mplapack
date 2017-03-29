/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaqp2.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlaqp2(INTEGER m, INTEGER n, INTEGER offset, REAL * A, INTEGER lda, INTEGER * jpvt, REAL * tau, REAL * vn1, REAL * vn2, REAL * work)
{
    INTEGER i, j, mn;
    REAL aii;
    INTEGER pvt;
    REAL temp, mtemp1;
    REAL temp2, tol3z;
    INTEGER offpi, itemp;
    REAL One = 1.0, Zero = 0.0;

    mn = min(m - offset, n);
    tol3z = sqrt(Rlamch("Epsilon"));
//Compute factorization.
    for (i = 0; i < mn; i++) {
	offpi = offset + i;
//Determine ith pivot column and swap if necessary.
	pvt = i - 1 + iRamax(n - i + 1, &vn1[i], 1);
	if (pvt != i) {
	    Rswap(m, &A[pvt * lda], 1, &A[i * lda], 1);
	    itemp = jpvt[pvt];
	    jpvt[pvt] = jpvt[i];
	    jpvt[i] = itemp;
	    vn1[pvt] = vn1[i];
	    vn2[pvt] = vn2[i];
	}
//Generate elementary reflector H(i).

	if (offpi < m) {
	    Rlarfg(m - offpi + 1, &A[offpi + i * lda], &A[offpi + 1 + i * lda], 1, &tau[i]);
	} else {
	    Rlarfg(1, &A[m + i * lda], &A[m + i * lda], 1, &tau[i]);
	}

	if (i < n) {
// Apply H(i)' to A(offset+i:m,i+1:n) from the left.
	    aii = A[offpi + i * lda];
	    A[offpi + i * lda] = One;
	    Rlarf("Left", m - offpi + 1, n - i, &A[offpi + i * lda], 1, tau[i], &A[offpi + (i + 1) * lda], lda, &work[0]);
	    A[offpi + i * lda] = aii;
	}
//Update partial column norms.
	for (j = i + 1; j < n; j++) {
	    if (vn1[j] != Zero) {

//NOTE: The following 4 lines follow from the analysis in
//Lapack Working Note 176.
		mtemp1 = abs(A[offpi + j * lda]) / vn1[j];
		temp = One - mtemp1 * mtemp1;
		temp = max(temp, Zero);
		temp2 = temp * ((vn1[j] / vn2[j]) * (vn1[j] / vn2[j]));
		if (temp2 <= tol3z) {
		    if (offpi < m) {
			vn1[j] = Rnrm2(m - offpi, &A[offpi + 1 + j * lda], 1);
			vn2[j] = vn1[j];
		    } else {
			vn1[j] = Zero;
			vn2[j] = Zero;
		    }
		} else {
		    vn1[j] *= sqrt(temp);
		}
	    }
	}
    }
    return;
}
