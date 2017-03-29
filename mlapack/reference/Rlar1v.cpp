/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlar1v.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MFALSE 0
#define MTRUE 1

void
Rlar1v(INTEGER n, INTEGER b1, INTEGER bn, REAL lambda, REAL * d, REAL * l,
       REAL * ld, REAL * lld, REAL pivmin, REAL gaptol,
       REAL * z, INTEGER wantnc, INTEGER * negcnt, REAL * ztz, REAL * mingma, INTEGER * r, INTEGER * isuppz, REAL * nrminv, REAL * resid, REAL * rqcorr, REAL * work)
{
    INTEGER i;
    REAL s;
    INTEGER r1, r2;
    REAL eps, tmp;
    INTEGER neg1, neg2, indp, inds;
    REAL dplus;
    INTEGER indlpl, indumn;
    REAL dminus;
    INTEGER sawnan1, sawnan2;
    REAL Zero = 0.0, One = 1.0;

    eps = Rlamch("Precision");
    if (*r == 0) {
	r1 = b1;
	r2 = bn;
    } else {
	r1 = *r;
	r2 = *r;
    }
//Storage for LPLUS
    indlpl = 0;
//Storage for UMINUS
    indumn = n;
    inds = (n * 2) + 1;
    indp = n * 3 + 1;
    if (b1 == 1) {
	work[inds] = Zero;
    } else {
	work[inds + b1 - 1] = lld[b1 - 1];
    }
//Compute the stationary transform (using the differential form)
//until the index R2
    sawnan1 = MFALSE;
    neg1 = 0;
    s = work[inds + b1 - 1] - lambda;
    for (i = b1; i <= r1 - 1; i++) {
	dplus = d[i] + s;
	work[indlpl + i] = ld[i] / dplus;
	if (dplus < Zero) {
	    ++neg1;
	}
	work[inds + i] = s * work[indlpl + i] * l[i];
	s = work[inds + i] - lambda;
    }
    sawnan1 = Risnan(s);
    if (sawnan1) {
	goto L60;
    }
    for (i = r1; i <= r2 - 1; i++) {
	dplus = d[i] + s;
	work[indlpl + i] = ld[i] / dplus;
	work[inds + i] = s * work[indlpl + i] * l[i];
	s = work[inds + i] - lambda;
    }
    sawnan1 = Risnan(s);
  L60:
    if (sawnan1) {
//Runs a slower version of the above loop if a NaN is detected
	neg1 = 0;
	s = work[inds + b1 - 1] - lambda;
	for (i = b1; i <= r1 - 1; i++) {
	    dplus = d[i] + s;
	    if (abs(dplus) < pivmin) {
		dplus = -(pivmin);
	    }
	    work[indlpl + i] = ld[i] / dplus;
	    if (dplus < Zero) {
		++neg1;
	    }
	    work[inds + i] = s * work[indlpl + i] * l[i];
	    if (work[indlpl + i] == Zero) {
		work[inds + i] = lld[i];
	    }
	    s = work[inds + i] - lambda;
	}
	for (i = r1; i <= r2 - 1; i++) {
	    dplus = d[i] + s;
	    if (abs(dplus) < pivmin) {
		dplus = -(pivmin);
	    }
	    work[indlpl + i] = ld[i] / dplus;
	    work[inds + i] = s * work[indlpl + i] * l[i];
	    if (work[indlpl + i] == Zero) {
		work[inds + i] = lld[i];
	    }
	    s = work[inds + i] - lambda;
	}
    }
//Compute the progressive transform (using the differential form)
//until the index R1
    sawnan2 = MFALSE;
    neg2 = 0;
    work[indp + bn - 1] = d[bn] - lambda;
    for (i = bn - 1; i >= r1; i--) {
	dminus = lld[i] + work[indp + i];
	tmp = d[i] / dminus;
	if (dminus < Zero) {
	    ++neg2;
	}
	work[indumn + i] = l[i] * tmp;
	work[indp + i - 1] = work[indp + i] * tmp - lambda;
    }
    tmp = work[indp + r1 - 1];
    sawnan2 = Risnan(tmp);
    if (sawnan2) {
//Runs a slower version of the above loop if a NaN is detected
	neg2 = 0;
	for (i = bn - 1; i >= r1; i--) {
	    dminus = lld[i] + work[indp + i];
	    if (abs(dminus) < pivmin) {
		dminus = -(pivmin);
	    }
	    tmp = d[i] / dminus;
	    if (dminus < Zero) {
		++neg2;
	    }
	    work[indumn + i] = l[i] * tmp;
	    work[indp + i - 1] = work[indp + i] * tmp - lambda;
	    if (tmp == Zero) {
		work[indp + i - 1] = d[i] - lambda;
	    }
	}
    }
//Find the index (from R1 to R2) of the largest (in magnitude)
//diagonal element of the inverse
    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (*mingma < Zero) {
	++neg1;
    }
    if (wantnc) {
	*negcnt = neg1 + neg2;
    } else {
	*negcnt = -1;
    }
    if (abs(*mingma) == Zero) {
	*mingma = eps * work[inds + r1 - 1];
    }
    *r = r1;
    for (i = r1; i <= r2 - 1; i++) {
	tmp = work[inds + i] + work[indp + i];
	if (tmp == Zero) {
	    tmp = eps * work[inds + i];
	}
	if (abs(tmp) <= abs(*mingma)) {
	    *mingma = tmp;
	    *r = i + 1;
	}
    }
//Compute the FP vector: solve N^T v = e_r
    isuppz[1] = b1;
    isuppz[2] = bn;
    z[*r] = One;
    *ztz = One;
//Compute the FP vector upwards from R
    if (!sawnan1 && !sawnan2) {
	for (i = *r - 1; i >= b1; i--) {
	    z[i] = -(work[indlpl + i] * z[i + 1]);
	    if ((abs(z[i]) + abs(z[i + 1])) * abs(ld[i]) < gaptol) {
		z[i] = Zero;
		isuppz[1] = i + 1;
		goto L220;
	    }
	    *ztz += z[i] * z[i];
	}
      L220:
	;
    } else {
//Run slower loop if NaN occurred.
	for (i = *r - 1; i >= b1; i--) {
	    if (z[i + 1] == Zero) {
		z[i] = -(ld[i + 1] / ld[i]) * z[i + 2];
	    } else {
		z[i] = -(work[indlpl + i] * z[i + 1]);
	    }
	    if ((abs(z[i]) + abs(z[i + 1])) * abs(ld[i]) < gaptol) {
		z[i] = Zero;
		isuppz[1] = i + 1;
		goto L240;
	    }
	    *ztz += z[i] * z[i];
	}
      L240:
	;
    }
//Compute the FP vector downwards from R in blocks of size BLKSIZ
    if (!sawnan1 && !sawnan2) {
	for (i = *r; i <= bn - 1; i++) {
	    z[i + 1] = -(work[indumn + i] * z[i]);
	    if ((abs(z[i]) + abs(z[i + 1])) * abs(ld[i]) < gaptol) {
		z[i + 1] = Zero;
		isuppz[2] = i;
		goto L260;
	    }
	    *ztz = *ztz + z[i + 1] * z[i + 1];

	}
      L260:
	;
    } else {
//Run slower loop if NaN occurred.
	for (i = *r; i <= bn - 1; i++) {
	    if (z[i] == Zero) {
		z[i + 1] = -(ld[i - 1] / ld[i]) * z[i - 1];
	    } else {
		z[i + 1] = -(work[indumn + i] * z[i]);
	    }
	    if ((abs(z[i]) + abs(z[i + 1])) * abs(ld[i]) < gaptol) {
		z[i + 1] = Zero;
		isuppz[2] = i;
		goto L280;
	    }
	    *ztz = *ztz + z[i + 1] * z[i + 1];
	}
      L280:
	;
    }
//Compute quantities for convergence test
    tmp = One / *ztz;
    *nrminv = sqrt(tmp);
    *resid = abs(*mingma) * (*nrminv);
    *rqcorr = *mingma * tmp;
    return;
}
