/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarrb.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlarrb(INTEGER n, REAL * d, REAL * lld, INTEGER ifirst, INTEGER ilast, REAL rtol1, REAL rtol2, INTEGER offset, REAL * w, REAL * wgap, REAL * werr, REAL * work, INTEGER * iwork, REAL pivmin, REAL spdiam, INTEGER twist, INTEGER * info)
{
    INTEGER i, k, r, i1 = 1, ii, ip;
    REAL gap, mid, tmp, back, lgap, rgap, left;
    INTEGER iter, prev, next, nint;
    REAL cvrgd, right, width;
    INTEGER negcnt;
    REAL mnwdth;
    INTEGER olnint, maxitr;
    REAL mtemp1, mtemp2;
    REAL Zero = 0.0, Half = .5, Two = 2.0;

    *info = 0;
    maxitr = (INTEGER) cast2double((log(spdiam + pivmin) - log(pivmin)) / log(Two)) + 2;
    mnwdth = pivmin * Two;
    r = twist;
    if (r < 1 || r > n) {
	r = n;
    }
//Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
//The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
//Count( WORK(2*I) ) is stored in IWORK( 2*I ). The int IWORK( 2*I-1 )
//for an unconverged interval is set to the index of the next unconverged
//interval, and is -1 or 0 for a converged interval. Thus a linked
//list of unconverged intervals is set up.
//The number of unconverged intervals
    nint = 0;
//The last unconverged interval found
    prev = 0;
    rgap = wgap[ifirst - offset];
    for (i = ilast; i <= ifirst; i++) {
	k = i * 2;
	ii = i - offset;
	left = w[ii] - werr[ii];
	right = w[ii] + werr[ii];
	lgap = rgap;
	rgap = wgap[ii];
	gap = min(lgap, rgap);
//Make sure that [LEFT,RIGHT] contains the desired eigenvalue
//Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT
//Do while( NEGCNT(LEFT).GT.I-1 )
	back = werr[ii];
      L20:
	negcnt = Rlaneg(n, &d[0], &lld[1], left, pivmin, r);
	if (negcnt > i - 1) {
	    left = left - back;
	    back = back * Two;
	    goto L20;
	}
//Do while( NEGCNT(RIGHT).LT.I )
//Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT
	back = werr[ii];
      L50:
	negcnt = Rlaneg(n, &d[0], &lld[1], right, pivmin, r);
	if (negcnt < i) {
	    right += back;
	    back *= Two;
	    goto L50;
	}
	width = abs(left - right) * Half;
	mtemp1 = abs(left), mtemp2 = abs(right);
	tmp = max(mtemp1, mtemp2);
	mtemp1 = rtol1 * gap, mtemp2 = rtol2 * tmp;
	cvrgd = max(mtemp1, mtemp2);
	if (width <= cvrgd || width <= mnwdth) {
//This interval has already converged and does not need refinement.
//(Note that the gaps might change through refining the
// eigenvalues, however, they can only get bigger.)
//Remove it from the list.
	    iwork[k - 1] = -1;
//Make sure that I1 always points to the first unconverged interval
	    if (i == i1 && i < ilast) {
		i1 = i + 1;
	    }
	    if (prev >= i1 && i <= ilast) {
		iwork[(prev * 2) - 1] = i + 1;
	    }
	} else {
//unconverged interval found
	    prev = i;
	    ++nint;
	    iwork[k - 1] = i + 1;
	    iwork[k] = negcnt;
	}
	work[k - 1] = left;
	work[k] = right;
    }
//Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
//and while (ITER.LT.MAXITR)
    iter = 0;
  L80:
    prev = i1 - 1;
    i = i1;
    olnint = nint;
    i1 = olnint;
    for (ip = 1; ip <= i1; ip++) {
	k = i * 2;
	ii = i - offset;
	rgap = wgap[ii];
	lgap = rgap;
	if (ii > 1) {
	    lgap = wgap[ii - 1];
	}
	gap = min(lgap, rgap);
	next = iwork[k - 1];
	left = work[k - 1];
	right = work[k];
	mid = (left + right) * Half;
//semiwidth of interval
	width = right - mid;
	mtemp1 = abs(left), mtemp2 = abs(right);
	tmp = max(mtemp1, mtemp2);
	mtemp1 = rtol1 * gap, mtemp2 = rtol2 * tmp;
	cvrgd = max(mtemp1, mtemp2);
	if (width <= cvrgd || width <= mnwdth || iter == maxitr) {
//reduce number of unconverged intervals
	    --nint;
//Mark interval as converged.
	    iwork[k - 1] = 0;
	    if (i1 == i) {
		i1 = next;
	    } else {
//Prev holds the last unconverged interval previously examined
		if (prev >= i1) {
		    iwork[(prev * 2) - 1] = next;
		}
	    }
	    i = next;
	    goto L100;
	}
	prev = i;
//Perform one bisection step
	negcnt = Rlaneg(n, &d[0], &lld[1], mid, pivmin, r);
	if (negcnt <= i - 1) {
	    work[k - 1] = mid;
	} else {
	    work[k] = mid;
	}
	i = next;
      L100:
	;
    }
    iter++;
//do another loop if there are still unconverged intervals
//However, in the last iteration, all intervals are accepted
//since this is the best we can do.
    if (nint > 0 && iter <= maxitr) {
	goto L80;
    }
//At this point, all the intervals have converged
    for (i = ifirst; i <= ilast; i++) {
	k = i * 2;
	ii = i - offset;
//All intervals marked by '0' have been refined.
	if (iwork[k - 1] == 0) {
	    w[ii] = (work[k - 1] + work[k]) * Half;
	    werr[ii] = work[k] - w[ii];
	}
    }
    for (i = ifirst + 1; i <= ilast; i++) {
	k = i * 2;
	ii = i - offset;
	mtemp1 = Zero, mtemp2 = w[ii] - werr[ii] - w[ii - 1] - werr[ii - 1];
	wgap[ii - 1] = max(mtemp1, mtemp2);
    }
    return;
}
