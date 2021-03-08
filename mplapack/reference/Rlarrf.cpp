/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarrf.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rlarrf(INTEGER n, REAL * d, REAL * l, REAL * ld, INTEGER clstrt, INTEGER clend, REAL * w,
	    REAL * wgap, REAL * werr, REAL spdiam, REAL clgapl, REAL clgapr, REAL pivmin, REAL * sigma, REAL * dplus, REAL * lplus, REAL * work, INTEGER * info)
{
    INTEGER i;
    REAL s, bestshift, smlgrowth, eps, tmp, max1, max2, rrr1, rrr2, znm2, growthbound, fail, fact, oldp;
    INTEGER indx = 0;
    REAL prod;
    INTEGER ktry;
    REAL fail2, avgap, ldmax, rdmax;
    INTEGER shift;
    INTEGER dorrr1;
    REAL ldelta;
    INTEGER nofail;
    REAL mingap, lsigma, rdelta;
    INTEGER forcer;
    REAL rsigma, clwdth;
    INTEGER sawnan1, sawnan2, tryrrr1;
    REAL Quarter = .25, One = 1.0, Two = 2.0, Four = 4.0, Eight = 8.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    fact = Two;
    eps = Rlamch("Precision");
    shift = 0;
    forcer = MFALSE;
//Note that we cannot guarantee that for any of the shifts tried,
//the factorization has a small or even moderate element growth.
//There could be Ritz values at both ends of the cluster and despite
//backing off, there are examples where all factorizations tried
//(in IEEE mode, allowing zero pivots & infinities) have INFINITE
//element growth.
//For this reason, we should use PIVMIN in this subroutine so that at
//least the L D L^T factorization exists. It can be checked afterwards
//whether the element growth caused bad residuals/orthogonality.
//Decide whether the code should accept the best among all
//representations despite large element growth or signal INFO=1
    nofail = MTRUE;

//Compute the average gap length of the cluster
    clwdth = abs(w[clend] - w[clstrt]) + werr[clend] + werr[clstrt];
    avgap = clwdth / (clend - clstrt);
    mingap = min(clgapl, clgapr);
//Initial values for shifts to both ends of cluster
    mtemp1 = w[clstrt], mtemp2 = w[clend];
    lsigma = min(mtemp1, mtemp2) - werr[clstrt];
    mtemp1 = w[clstrt], mtemp2 = w[clend];
    rsigma = max(mtemp1, mtemp2) + werr[clend];
//Use a small fudge to make sure that we really shift to the outside
    lsigma = lsigma - abs(lsigma) * Four * eps;
    rsigma = rsigma + abs(rsigma) * Four * eps;
//Compute upper bounds for how much to back off the initial shifts
    ldmax = mingap * Quarter + pivmin * Two;
    rdmax = mingap * Quarter + pivmin * Two;
    mtemp1 = avgap, mtemp2 = wgap[clstrt];
    ldelta = max(mtemp1, mtemp2) / fact;
    mtemp1 = avgap, mtemp2 = wgap[clend - 1];
    rdelta = max(mtemp1, mtemp2) / fact;
//Initialize the record of the best representation found
    s = Rlamch("S");
    smlgrowth = One / s;
    fail = (n - 1) * mingap / (spdiam * eps);
    fail2 = (n - 1) * mingap / (spdiam * sqrt(eps));
    bestshift = lsigma;
//while (KTRY <= KTRYMAX)
    ktry = 0;
    growthbound = spdiam * Eight;
  L5:
    sawnan1 = MFALSE;
    sawnan2 = MFALSE;
//Ensure that we do not back off too much of the initial shifts
    ldelta = min(ldmax, ldelta);
    rdelta = min(rdmax, rdelta);
//Compute the element growth when shifting to both ends of the cluster
//accept the shift if there is no element growth at one of the two ends
//Left end
    s = -lsigma;
    dplus[1] = d[1] + s;
    if (abs(dplus[1]) < pivmin) {
	dplus[1] = -(pivmin);
//Need to set SAWNAN1 because refined RRR test should not be used
//in this case
	sawnan1 = MTRUE;
    }
    max1 = abs(dplus[1]);
    for (i = 0; i < n - 1; i++) {
	lplus[i] = ld[i] / dplus[i];
	s = s * lplus[i] * l[i] - lsigma;
	dplus[i + 1] = d[i + 1] + s;
	if (abs(dplus[i + 1]) < pivmin) {
	    dplus[i + 1] = -(pivmin);
//Need to set SAWNAN1 because refined RRR test should not be used
//in this case
	    sawnan1 = MTRUE;
	}
	mtemp1 = max1, mtemp2 = abs(dplus[i + 1]);
	max1 = max(mtemp1, mtemp2);
    }
    sawnan1 = sawnan1 || Risnan(max1);
    if (forcer || (max1 <= growthbound && !sawnan1)) {
	*sigma = lsigma;
	shift = 1;
	goto L100;
    }
//Right end
    s = -rsigma;
    work[1] = d[1] + s;
    if (abs(work[1]) < pivmin) {
	work[1] = -(pivmin);
//Need to set SAWNAN2 because refined RRR test should not be used
//in this case
	sawnan2 = MTRUE;
    }
    max2 = abs(work[1]);
    for (i = 0; i < n - 1; i++) {
	work[n + i] = ld[i] / work[i];
	s = s * work[n + i] * l[i] - rsigma;
	work[i + 1] = d[i + 1] + s;
	if (abs(work[i + 1]) < pivmin) {
	    work[i + 1] = -(pivmin);
//Need to set SAWNAN2 because refined RRR test should not be used
//in this case
	    sawnan2 = MTRUE;
	}
	mtemp1 = max2, mtemp2 = abs(work[i + 1]);
	max2 = max(mtemp1, mtemp2);
    }
    sawnan2 = sawnan2 || Risnan(max2);
    if (forcer || (max2 <= growthbound && !sawnan2)) {
	*sigma = rsigma;
	shift = 2;
	goto L100;
    }
//If we are at this point, both shifts led to too much element growth
//Record the better of the two shifts (provided it didn't lead to NaN)
    if (sawnan1 && sawnan2) {
//both MAX1 and MAX2 are NaN
	goto L50;
    } else {
	if (!sawnan1) {
	    indx = 1;
	    if (max1 <= smlgrowth) {
		smlgrowth = max1;
		bestshift = lsigma;
	    }
	}
	if (!sawnan2) {
	    if (sawnan1 || max2 <= max1) {
		indx = 2;
	    }
	    if (max2 <= smlgrowth) {
		smlgrowth = max2;
		bestshift = rsigma;
	    }
	}
    }
//If we are here, both the left and the right shift led to
//element growth. If the element growth is moderate, then
//we may still accept the representation, if it passes a
//refined test for RRR. This test supposes that no NaN occurred.
//Moreover, we use the refined RRR test only for isolated clusters.
    if (clwdth < mingap / 128 && min(max1, max2) < fail2 && !sawnan1 && !sawnan2) {
	dorrr1 = MTRUE;
    } else {
	dorrr1 = MFALSE;
    }
    tryrrr1 = MTRUE;
    if (tryrrr1 && dorrr1) {
	if (indx == 1) {
	    tmp = abs(dplus[n]);
	    znm2 = One;
	    prod = One;
	    oldp = One;
	    for (i = n - 1; i >= 1; i--) {
		if (prod <= eps) {
		    prod = dplus[i + 1] * work[n + i + 1] / (dplus[i] * work[n + i]) * oldp;
		} else {
		    prod = prod * abs(work[n + i]);
		}
		oldp = prod;
		znm2 = znm2 + prod * prod;
		mtemp1 = tmp, mtemp2 = abs(dplus[i] * prod);
		tmp = max(mtemp1, mtemp2);
	    }
	    rrr1 = tmp / (spdiam * sqrt(znm2));
	    if (rrr1 <= Eight) {
		*sigma = lsigma;
		shift = 1;
		goto L100;
	    }
	} else if (indx == 2) {
	    tmp = abs(work[n]);
	    znm2 = One;
	    prod = One;
	    oldp = One;
	    for (i = n - 1; i >= 1; i--) {
		if (prod <= eps) {
		    prod = work[i + 1] * lplus[i + 1] / (work[i] * lplus[i]) * oldp;
		} else {
		    prod *= abs(lplus[i]);
		}
		oldp = prod;
		znm2 = znm2 + prod * prod;
		mtemp1 = tmp, mtemp2 = abs(work[i] * prod);
		tmp = max(mtemp1, mtemp2);
	    }
	    rrr2 = tmp / (spdiam * sqrt(znm2));
	    if (rrr2 <= Eight) {
		*sigma = rsigma;
		shift = 2;
		goto L100;
	    }
	}
    }
  L50:
    if (ktry < 1) {
//If we are here, both shifts failed also the RRR test.
//Back off to the outside
	mtemp1 = lsigma - ldelta, mtemp2 = lsigma - ldmax;
	lsigma = max(mtemp1, mtemp2);
	mtemp1 = rsigma + rdelta, mtemp2 = rsigma + rdmax;
	rsigma = min(mtemp1, mtemp2);
	ldelta = ldelta * Two;
	rdelta = rdelta * Two;
	ktry++;
	goto L5;
    } else {
//None of the representations investigated satisfied our
//criteria. Take the best one we found.
	if (smlgrowth < fail || nofail) {
	    lsigma = bestshift;
	    rsigma = bestshift;
	    forcer = MTRUE;
	    goto L5;
	} else {
	    *info = 1;
	    return;
	}
    }
  L100:
    if (shift == 1) {
    } else if (shift == 2) {
//store new L and D back into DPLUS, LPLUS
	Rcopy(n, &work[0], 1, &dplus[1], 1);
	Rcopy(n - 1, &work[n + 1], 1, &lplus[1], 1);
    }
    return;
}
