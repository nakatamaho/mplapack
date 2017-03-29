/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarft.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clarft(const char *direct, const char *storev, INTEGER n, INTEGER k, COMPLEX * v, INTEGER ldv, COMPLEX * tau, COMPLEX * t, INTEGER ldt)
{
    INTEGER i, j;
    COMPLEX vii;
    REAL Zero = 0.0, One = 1.0;

    //Quick return if possible
    if (n == 0) {
	return;
    }
    if (Mlsame(direct, "F")) {
	for (i = 1; i <= k; i++) {
	    if (tau[i - 1] == Zero) {
		//H(i)  =  I
		for (j = 1; j <= i; j++) {
		    t[(j - 1) + (i - 1) * ldt] = Zero;
		}
	    } else {
		//general case
		vii = v[(i - 1) + (i - 1) * ldv];
		v[(i - 1) + (i - 1) * ldv] = One;
		if (Mlsame(storev, "C")) {
		    // T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
		    Cgemv("Conjugate Transpose", n - i + 1, i - 1, -tau[i - 1], &v[(i - 1) + 0 * ldv], ldv, &v[(i - 1) + (i - 1) * ldv], 1, Zero, &t[0 + (i - 1) * ldt], 1);
		} else {
		    //T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
		    if (i < n) {
			Clacgv(n - i, &v[(i - 1) + i * ldv], ldv);
		    }
		    Cgemv("No transpose", i - 1, n - i + 1, -tau[i - 1], &v[0 + (i - 1) * ldv], ldv, &v[(i - 1) + (i - 1) * ldv], ldv, Zero, &t[0 + (i - 1) * ldt], 1);
		    if (i < n) {
			Clacgv(n - i, &v[(i - 1) + i * ldv], ldv);
		    }
		}
		v[(i - 1) + (i - 1) * ldv] = vii;
		//T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
		Ctrmv("Upper", "No transpose", "Non-unit", i - 1, t, ldt, &t[0 + (i - 1) * ldt], 1);
		t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
	    }
	}
    } else {
	for (i = k; i >= 1; i--) {
	    if (tau[i - 1] == Zero) {
		//H(i)  =  I
		for (j = i; j < k; j++) {
		    t[(j - 1) + (i - 1) * ldt] = Zero;
		}
	    } else {
		//general case
		if (i < k) {
		    if (Mlsame(storev, "C")) {
			vii = v[(n - k + i - 1) + (i - 1) * ldv];
			v[(n - k + i - 1) + (i - 1) * ldv] = One;
			//T(i+1:k,i) := - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
			Cgemv("Conjugate Transpose", n - k + i, k - i, -tau[i - 1], &v[0 + i * ldv], ldv, &v[0 + (i - 1) * ldv], 1, Zero, &t[i + (i - 1) * ldt], 1);
			v[(n - k + i - 1) + (i - 1) * ldv] = vii;
		    } else {
			vii = v[(i - 1) + (n - k + i - 1) * ldv];
			v[(i - 1) + (n - k + i - 1) * ldv] = One;
			//T(i+1:k,i) := - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
			Clacgv(n - k + i - 1, &v[i - 1], ldv);
			Cgemv("No transpose", k - i, n - k + i, -tau[i - 1], &v[i + 0 * ldv], ldv, &v[(i - 1) + 0 * ldv], ldv, Zero, &t[i + (i - 1) * ldt], 1);
			Clacgv(n - k + i - 1, &v[i - 1], ldv);
			v[(i - 1) + (n - k + i - 1) * ldv] = vii;
		    }
		    //T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
		    Ctrmv("Lower", "No transpose", "Non-unit", k - i, &t[i + i * ldt], ldt, &t[i + (i - 1) * ldt], 1);
		}
		t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
	    }
	}
    }
    return;
}
