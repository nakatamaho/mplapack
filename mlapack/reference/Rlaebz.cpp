/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaebz.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rlaebz(INTEGER ijob, INTEGER nitmax, INTEGER n, INTEGER mmax, INTEGER minp, INTEGER nbmin,
       REAL abstol, REAL reltol, REAL pivmin, REAL * d,
       REAL * e, REAL * e2, INTEGER * nval, REAL * AB, REAL * c, INTEGER * mout, INTEGER * nab, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER j, kf, ji, kl, jp, jit;
    REAL tmp1, tmp2;
    INTEGER itmp1, itmp2, kfnew, klnew;
    REAL Half = 0.5, Zero = 0.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    if (ijob < 1 || ijob > 3) {
	*info = -1;
	return;
    }
//Initialize NAB
    if (ijob == 1) {
//Compute the number of eigenvalues in the initial intervals.
	mout = 0;
	for (ji = 0; ji <= minp; ji++) {
	    for (jp = 1; jp <= 2; jp++) {
		tmp1 = d[1] - AB[ji + jp * mmax];
		if (abs(tmp1) < pivmin) {
		    tmp1 = -(pivmin);
		}
		nab[ji + jp * mmax] = 0;
		if (tmp1 <= Zero) {
		    nab[ji + jp * mmax] = 1;
		}
		for (j = 2; j <= n; j++) {
		    tmp1 = d[j] - e2[j - 1] / tmp1 - AB[ji + jp * mmax];
		    if (abs(tmp1) < pivmin) {
			tmp1 = -(pivmin);
		    }
		    if (tmp1 <= Zero) {
			++nab[ji + jp * mmax];
		    }
		}
	    }
	    mout = mout + nab[ji + (mmax << 1)] - nab[ji + mmax];
	}
	return;
    }
//Initialize for loop
//KF and KL have the following meaning:
//   Intervals 1,...,KF-1 have converged.
//   Intervals KF,...,KL  still need to be refined.
    kf = 1;
    kl = minp;
//If IJOB=2, initialize C.
//If IJOB=3, use the user-supplied starting point.
    if (ijob == 2) {
	for (ji = 0; ji <= minp; ji++) {
	    c[ji] = (AB[ji + mmax] + AB[ji + (mmax << 1)]) * Half;
	}
    }
//Iteration loop
    for (jit = 1; jit <= nitmax; jit++) {
//Loop over intervals
	if (kl - kf + 1 >= nbmin && nbmin > 0) {
//Begin of Parallel Version of the loop
	    for (ji = kf; ji <= kl; ji++) {
//Compute N(c), the number of eigenvalues less than c
		work[ji] = d[1] - c[ji];
		iwork[ji] = 0;
		if (work[ji] <= pivmin) {
		    iwork[ji] = 1;
		    mtemp1 = work[ji], mtemp2 = -(pivmin);
		    work[ji] = min(mtemp1, mtemp2);
		}
		for (j = 2; j <= n; j++) {
		    work[ji] = d[j] - e2[j - 1] / work[ji] - c[ji];
		    if (work[ji] <= pivmin) {
			iwork[ji]++;
			mtemp1 = work[ji], mtemp2 = -(pivmin);
			work[ji] = min(mtemp1, mtemp2);
		    }
		}
	    }
	    if (ijob <= 2) {
//IJOB=2: Choose all intervals containing eigenvalues.
		klnew = kl;
		for (ji = kf; ji <= kl; ji++) {
//Insure that N(w) is monotone
		    iwork[ji] = min(nab[ji + (mmax << 1)], max(nab[ji + mmax], iwork[ji]));
//Update the Queue -- add intervals if both halves
//contain eigenvalues.
		    if (iwork[ji] == nab[ji + (mmax << 1)]) {
//No eigenvalue in the upper interval:
//just use the lower interval.
			AB[ji + (mmax << 1)] = c[ji];
		    } else if (iwork[ji] == nab[ji + mmax]) {
//No eigenvalue in the lower interval:
//just use the upper interval.
			AB[ji + mmax] = c[ji];
		    } else {
			klnew++;
			if (klnew <= mmax) {
//Eigenvalue in both intervals -- add upper to
//queue.
			    AB[klnew + (mmax << 1)] = AB[ji + (mmax << 1)];
			    nab[klnew + (mmax << 1)] = nab[ji + (mmax << 1)];
			    AB[klnew + mmax] = c[ji];
			    nab[klnew + mmax] = iwork[ji];
			    AB[ji + (mmax << 1)] = c[ji];
			    nab[ji + (mmax << 1)] = iwork[ji];
			} else {
			    *info = mmax + 1;
			}
		    }
		}
		if (*info != 0) {
		    return;
		}
		kl = klnew;
	    } else {
//IJOB=3: Binary search.  Keep only the interval containing
//        w   s.t. N(w) = NVAL
		for (ji = kf; ji <= kl; ji++) {
		    if (iwork[ji] <= nval[ji]) {
			AB[ji + mmax] = c[ji];
			nab[ji + mmax] = iwork[ji];
		    }
		    if (iwork[ji] >= nval[ji]) {
			AB[ji + (mmax << 1)] = c[ji];
			nab[ji + (mmax << 1)] = iwork[ji];
		    }

		}
	    }
	} else {
//End of Parallel Version of the loop
//Begin of Serial Version of the loop
	    klnew = kl;
	    for (ji = kf; ji <= kl; ji++) {
//Compute N(w), the number of eigenvalues less than w
		tmp1 = c[ji];
		tmp2 = d[1] - tmp1;
		itmp1 = 0;
		if (tmp2 <= pivmin) {
		    itmp1 = 1;
		    mtemp1 = tmp2, mtemp2 = -(pivmin);
		    tmp2 = min(mtemp1, mtemp2);
		}
		for (j = 2; j <= n; j++) {
		    tmp2 = d[j] - e2[j - 1] / tmp2 - tmp1;
		    if (tmp2 <= pivmin) {
			itmp1++;
			mtemp1 = tmp2, mtemp2 = -(pivmin);
			tmp2 = min(mtemp1, mtemp2);
		    }
		}
		if (ijob <= 2) {
//IJOB=2: Choose all intervals containing eigenvalues.
//Insure that N(w) is monotone
		    itmp1 = min(nab[ji + (mmax << 1)], max(nab[ji + mmax], itmp1));
//Update the Queue -- add intervals if both halves
//contain eigenvalues.
		    if (itmp1 == nab[ji + (mmax << 1)]) {
//No eigenvalue in the upper interval:
//just use the lower interval.
			AB[ji + (mmax << 1)] = tmp1;
		    } else if (itmp1 == nab[ji + mmax]) {
//No eigenvalue in the lower interval:
//just use the upper interval.
			AB[ji + mmax] = tmp1;
		    } else if (klnew < mmax) {
//Eigenvalue in both intervals -- add upper to queue.
			klnew++;
			AB[klnew + (mmax << 1)] = AB[ji + (mmax << 1)];
			nab[klnew + (mmax << 1)] = nab[ji + (mmax << 1)];
			AB[klnew + mmax] = tmp1;
			nab[klnew + mmax] = itmp1;
			AB[ji + (mmax << 1)] = tmp1;
			nab[ji + (mmax << 1)] = itmp1;
		    } else {
			*info = mmax + 1;
			return;
		    }
		} else {
//IJOB=3: Binary search.  Keep only the interval
//        containing  w  s.t. N(w) = NVAL
		    if (itmp1 <= nval[ji]) {
			AB[ji + mmax] = tmp1;
			nab[ji + mmax] = itmp1;
		    }
		    if (itmp1 >= nval[ji]) {
			AB[ji + (mmax << 1)] = tmp1;
			nab[ji + (mmax << 1)] = itmp1;
		    }
		}
	    }
	    kl = klnew;
//End of Serial Version of the loop
	}
//Check for convergence
	kfnew = kf;
	for (ji = kf; ji <= kl; ji++) {
	    tmp1 = abs(AB[ji + 2 * mmax] - AB[ji + mmax]);
	    tmp2 = max(abs(AB[ji + 2 * mmax]), abs(AB[ji + mmax]));
	    mtemp1 = max(abstol, pivmin), mtemp2 = reltol * tmp2;
	    if (tmp1 < max(mtemp1, mtemp2)
		|| nab[ji + mmax] >= nab[ji + (mmax << 1)]) {
//Converged -- Swap with position KFNEW,
//             then increment KFNEW
		if (ji > kfnew) {
		    tmp1 = AB[ji + mmax];
		    tmp2 = AB[ji + (mmax << 1)];
		    itmp1 = nab[ji + mmax];
		    itmp2 = nab[ji + (mmax << 1)];
		    AB[ji + mmax] = AB[kfnew + mmax];
		    AB[ji + (mmax << 1)] = AB[kfnew + (mmax << 1)];
		    nab[ji + mmax] = nab[kfnew + mmax];
		    nab[ji + (mmax << 1)] = nab[kfnew + (mmax << 1)];
		    AB[kfnew + mmax] = tmp1;
		    AB[kfnew + (mmax << 1)] = tmp2;
		    nab[kfnew + mmax] = itmp1;
		    nab[kfnew + (mmax << 1)] = itmp2;
		    if (ijob == 3) {
			itmp1 = nval[ji];
			nval[ji] = nval[kfnew];
			nval[kfnew] = itmp1;
		    }
		}
		kfnew++;
	    }
	}
	kf = kfnew;
//Choose Midpoints
	for (ji = kf; ji <= kl; ji++) {
	    c[ji] = (AB[ji + mmax] + AB[ji + (mmax << 1)]) * Half;
	}
//If no more intervals to refine, quit.
	if (kf > kl) {
	    goto L140;
	}
    }
//Converged
  L140:
    *info = max(kl + 1 - kf, (INTEGER) 0);
    *mout = kl;
    return;
}
