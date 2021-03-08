/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cunm2r.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cunm2r(const char *side, const char *trans, INTEGER m, INTEGER n, INTEGER k, COMPLEX * A, INTEGER lda, COMPLEX * tau, COMPLEX * c, INTEGER ldc, COMPLEX * work, INTEGER * info)
{
    INTEGER i, i1, i2, i3, ic, jc, mi = 0, ni = 0, nq;
    COMPLEX aii;
    INTEGER left;
    COMPLEX taui;
    INTEGER notran;
    REAL One = 1.0;

    *info = 0;
    left = Mlsame(side, "L");
    notran = Mlsame(trans, "N");
//NQ is the order of Q
    if (left) {
	nq = m;
    } else {
	nq = n;
    }
    if (!left && !Mlsame(side, "R")) {
	*info = -1;
    } else if (!notran && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (k < 0 || k > nq) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, nq)) {
	*info = -7;
    } else if (ldc < max((INTEGER) 1, m)) {
	*info = -10;
    }
    if (*info != 0) {
	Mxerbla("Cunm2r", -(*info));
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0 || k == 0) {
	return;
    }
    if ((left && !notran) || (!left && notran)) {
	i1 = 1;
	i2 = k;
	i3 = 1;
    } else {
	i1 = k;
	i2 = 1;
	i3 = -1;
    }
    if (left) {
	ni = n;
	jc = 1;
    } else {
	mi = m;
	ic = 1;
    }
    for (i = i1; i <= i2; i += i3) {
	if (left) {
//H(i) or H(i)' is applied to C(i:m,1:n)
	    mi = m - i + 1;
	    ic = i;
	} else {
//H(i) or H(i)' is applied to C(1:m,i:n)
	    ni = n - i + 1;
	    jc = i;
	}
//Apply H(i) or H(i)'
	if (notran) {
	    taui = tau[i];
	} else {
	    taui = conj(tau[i]);
	}
	aii = A[i + i * lda];
	A[i + i * lda] = One;
	Clarf(side, mi, ni, &A[i + i * lda], 1, taui, &c[ic + jc * ldc], ldc, &work[0]);
	A[i + i * lda] = aii;
    }
    return;
}
