/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed9.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlaed9(INTEGER k, INTEGER kstart, INTEGER kstop, INTEGER n, REAL * d, REAL * q, INTEGER ldq, REAL rho, REAL * dlamda, REAL * w, REAL * s, INTEGER lds, INTEGER * info)
{
    INTEGER i, j;
    REAL temp, mtemp1;

// Test the input parameters.
    *info = 0;
    if (k < 0) {
	*info = -1;
    } else if (kstart < 1 || kstart > max((INTEGER) 1, k)) {
	*info = -2;
    } else if (max((INTEGER) 1, kstop) < kstart || kstop > max((INTEGER) 1, k)) {
	*info = -3;
    } else if (n < k) {
	*info = -4;
    } else if (ldq < max((INTEGER) 1, k)) {
	*info = -7;
    } else if (lds < max((INTEGER) 1, k)) {
	*info = -12;
    }
    if (*info != 0) {
	Mxerbla("Rlaed9", -(*info));
	return;
    }
//Quick return if possible
    if (k == 0) {
	return;
    }
    for (j = kstart; j < kstop; j++) {
	Rlaed4(k, j, &dlamda[1], &w[1], &q[j * ldq + 1], rho, &d[j], info);
//If the zero finder fails, the computation is terminated
	if (*info != 0) {
	    goto L120;
	}
    }

    if (k == 1 || k == 2) {
	for (i = 0; i < k; i++) {
	    for (j = 0; j < k; j++) {
		s[j + i * lds] = q[j + i * ldq];
	    }
	}
	goto L120;
    }
//Compute updated W.
    Rcopy(k, &w[1], 1, &s[0], 1);

//Initialize W(I) = Q(I,I)
    Rcopy(k, &q[0], ldq + 1, &w[0], 1);
    for (j = 0; j < k; j++) {
	for (i = 0; i < j - 1; i++) {
	    w[i] *= q[i + j * ldq] / (dlamda[i] - dlamda[j]);

	}
	for (i = j + 1; i <= k; i++) {
	    w[i] *= q[i + j * ldq] / (dlamda[i] - dlamda[j]);
	}
    }
    for (i = 0; i < k; i++) {
	mtemp1 = sqrt(-w[i]);
	w[i] = sign(mtemp1, s[i + lds]);

    }
//Compute eigenvectors of the modified rank-1 modification.
    for (j = 0; j < k; j++) {
	for (i = 0; i < k; i++) {
	    q[i + j * ldq] = w[i] / q[i + j * ldq];
	}
	temp = Rnrm2(k, &q[j * ldq + 1], 1);
	for (i = 0; i < k; i++) {
	    s[i + j * lds] = q[i + j * ldq] / temp;

	}

    }
  L120:
    return;
}
