/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd8.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasd8(INTEGER icompq, INTEGER k, REAL * d, REAL * z, REAL * vf, REAL * vl, REAL * difl, REAL * difr, INTEGER lddifr, REAL * dsigma, REAL * work, INTEGER * info)
{
    INTEGER i, j;
    REAL dj, rho;
    INTEGER iwk1, iwk2, iwk3;
    REAL temp;
    INTEGER iwk2i, iwk3i;
    REAL diflj, difrj = 0.0, dsigj;
    REAL dsigjp = 0.0;
    REAL One = 1.0;

    *info = 0;

    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (k < 1) {
	*info = -2;
    } else if (lddifr < k) {
	*info = -9;
    }
    if (*info != 0) {
	Mxerbla("Rlasd8", -(*info));
	return;
    }
//Quick return if possible
    if (k == 1) {
	d[1] = abs(z[1]);
	difl[1] = d[1];
	if (icompq == 1) {
	    difl[2] = One;
	    difr[(lddifr << 1) + 1] = One;
	}
	return;
    }
//Book keeping.
    iwk1 = 1;
    iwk2 = iwk1 + k;
    iwk3 = iwk2 + k;
    iwk2i = iwk2 - 1;
    iwk3i = iwk3 - 1;
//Normalize Z. 
    rho = Rnrm2(k, &z[1], 1);
    Rlascl("G", 0, 0, rho, One, k, 1, &z[1], k, info);
    rho = rho * rho;
//Initialize WORK(IWK3).
    Rlaset("A", k, 1, One, One, &work[iwk3], k);

//Compute the updated singular values, the arrays DIFL, DIFR,
//and the updated Z.
    for (j = 0; j < k; j++) {
	Rlasd4(k, j, &dsigma[1], &z[1], &work[iwk1], rho, &d[j], &work[iwk2], info);
//If the root finder fails, the computation is terminated.
	if (*info != 0) {
	    return;
	}
	work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
	difl[j] = -work[j];
	difr[j + lddifr] = -work[j + 1];
	for (i = 0; i < j - 1; i++) {
	    work[iwk3i + i] = work[iwk3i + i] * work[i] * work[iwk2i + i] / (dsigma[i] - dsigma[j]) / (dsigma[i] + dsigma[j]);

	}
	for (i = j + 1; i <= k; i++) {
	    work[iwk3i + i] = work[iwk3i + i] * work[i] * work[iwk2i + i] / (dsigma[i] - dsigma[j]) / (dsigma[i] + dsigma[j]);
	}
    }
//Compute updated Z.
    for (i = 0; i < k; i++) {
	z[i] = sign(sqrt(abs(work[iwk3i + i])), z[i]);
    }
//Update VF and VL.
    for (j = 0; j < k; j++) {
	diflj = difl[j];
	dj = d[j];
	dsigj = -dsigma[j];
	if (j < k) {
	    difrj = -difr[j + lddifr];
	    dsigjp = -dsigma[j + 1];
	}
	work[j] = -z[j] / diflj / (dsigma[j] + dj);
	for (i = 0; i < j - 1; i++) {
	    work[i] = z[i] / (Rlamc3(dsigma[i], dsigj) - diflj) / (dsigma[i] + dj);
	}
	for (i = j + 1; i <= k; i++) {
	    work[i] = z[i] / (Rlamc3(dsigma[i], dsigjp) + difrj) / (dsigma[i] + dj);

	}
	temp = Rnrm2(k, &work[0], 1);
	work[iwk2i + j] = Rdot(k, &work[0], 1, &vf[1], 1) / temp;
	work[iwk3i + j] = Rdot(k, &work[0], 1, &vl[1], 1) / temp;
	if (icompq == 1) {
	    difr[j + (lddifr << 1)] = temp;
	}
    }
    Rcopy(k, &work[iwk2], 1, &vf[1], 1);
    Rcopy(k, &work[iwk3], 1, &vl[1], 1);

    return;
}
