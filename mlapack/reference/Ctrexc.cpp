/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctrexc.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctrexc(const char *compq, INTEGER n, COMPLEX * t, INTEGER ldt, COMPLEX * q, INTEGER ldq, INTEGER ifst, INTEGER ilst, INTEGER * info)
{
    INTEGER k, m1, m2, m3;
    REAL cs;
    COMPLEX t11, t22, sn, temp;
    INTEGER wantq;

//Decode and test the input parameters.
    *info = 0;
    wantq = Mlsame(compq, "V");
    if (!Mlsame(compq, "N") && !wantq) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (ldt < max((INTEGER) 1, n)) {
	*info = -4;
    } else if (ldq < 1 || (wantq && ldq < max((INTEGER) 1, n))) {
	*info = -6;
    } else if (ifst < 1 || ifst > n) {
	*info = -7;
    } else if (ilst < 1 || ilst > n) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Ctrexc", -(*info));
	return;
    }
//Quick return if possible
    if (n == 1 || ifst == ilst) {
	return;
    }
    if (ifst < ilst) {
//Move the IFST-th diagonal element forward down the diagonal.
	m1 = 0;
	m2 = -1;
	m3 = 1;
    } else {
//Move the IFST-th diagonal element backward up the diagonal.
	m1 = -1;
	m2 = 0;
	m3 = -1;
    }
    for (k = ifst + m1; k <= ilst + m2; k = k + m3) {
//Interchange the k-th and (k+1)-th diagonal elements.
	t11 = t[k + k * ldt];
	t22 = t[k + 1 + (k + 1) * ldt];
//Determine the transformation to perform the interchange.
	Clartg(t[k + (k + 1) * ldt], t22 - t11, &cs, &sn, &temp);
//Apply transformation to the matrix T.
	if (k + 2 <= n) {
	    Crot(n - k - 1, &t[k + (k + 2) * ldt], ldt, &t[k + 1 + (k + 2) * ldt], ldt, cs, sn);
	}
	Crot(k - 1, &t[k * ldt + 1], 1, &t[(k + 1) * ldt + 1], 1, cs, conj(sn));
	t[k + k * ldt] = t22;
	t[k + 1 + (k + 1) * ldt] = t11;
	if (wantq) {
//Accumulate transformation in the matrix Q.
	    Crot(n, &q[k * ldq + 1], 1, &q[(k + 1) * ldq + 1], 1, cs, conj(sn));
	}
    }
    return;
}
