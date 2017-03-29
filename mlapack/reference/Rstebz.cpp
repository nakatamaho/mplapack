/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rstebz.cpp,v 1.6 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rstebz(const char *range, const char *order, INTEGER n, REAL vl, REAL vu, INTEGER il, INTEGER iu, REAL abstol, REAL * d,
	    REAL * e, INTEGER * m, INTEGER * nsplit, REAL * w, INTEGER * iblock, INTEGER * isplit, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER j, ib, jb, ie, je, nb;
    REAL gl;
    INTEGER im, in;
    REAL gu;
    INTEGER iw;
    REAL wl, wu;
    INTEGER nwl;
    REAL ulp, wlu = 0, wul = 0;
    INTEGER nwu;
    REAL tmp1, tmp2;
    INTEGER iend, ioff, iout, itmp1, jdisc;
    INTEGER iinfo;
    REAL atoli;
    INTEGER iwoff;
    REAL bnorm;
    INTEGER itmax;
    REAL wkill, rtoli, tnorm;
    INTEGER ibegin;
    INTEGER irange, idiscl;
    REAL safemn;
    INTEGER idumma;
    INTEGER idiscu, iorder;
    INTEGER ncnvrg;
    REAL pivmin;
    INTEGER toofew;
    REAL Two = 2.0, One = 1.0, Zero = 0.0, Fudge = 2.1, Half = 0.5;
    REAL mtemp1, mtemp2;

//Decode RANGE
    if (Mlsame(range, "A")) {
	irange = 1;
    } else if (Mlsame(range, "V")) {
	irange = 2;
    } else if (Mlsame(range, "I")) {
	irange = 3;
    } else {
	irange = 0;
    }
//Decode ORDER
    if (Mlsame(order, "B")) {
	iorder = 2;
    } else if (Mlsame(order, "E")) {
	iorder = 1;
    } else {
	iorder = 0;
    }
//Check for Errors
    if (irange <= 0) {
	*info = -1;
    } else if (iorder <= 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (irange == 2) {
	if (vl >= vu) {
	    *info = -5;
	}
    } else if (irange == 3 && (il < 1 || il > max((INTEGER) 1, n))) {
	*info = -6;
    } else if (irange == 3 && (iu < min(n, il) || iu > n)) {
	*info = -7;
    }
    if (*info != 0) {
	Mxerbla("Rstebz", -(*info));
	return;
    }
//Initialize error flags
    *info = 0;
    ncnvrg = MFALSE;
    toofew = MFALSE;
//Quick return if possible
    m = 0;
    if (n == 0) {
	return;
    }
//Simplifications:
    if (irange == 3 && il == 1 && iu == n) {
	irange = 1;
    }
//Get machine constants
//NB is the minimum vector length for vector bisection, or 0
//if only scalar is to be done.
    safemn = Rlamch("S");
    ulp = Rlamch("P");
    rtoli = ulp * Two;
    nb = iMlaenv(1, "Rstebz", " ", n, -1, -1, -1);
    if (nb <= 1) {
	nb = 0;
    }
//Special Case when N=1
    if (n == 1) {
	*nsplit = 1;
	isplit[1] = 1;
	if (irange == 2 && (vl >= d[1] || vu < d[1])) {
	    *m = 0;
	} else {
	    w[1] = d[1];
	    iblock[1] = 1;
	    *m = 1;
	}
	return;
    }
//Compute Splitting PoINTEGERs
    *nsplit = 1;
    work[n] = Zero;
    pivmin = One;
    for (j = 2; j <= n; j++) {
	tmp1 = e[j - 1] * e[j - 1];
	if (abs(d[j] * d[j - 1]) * (ulp * ulp) + safemn > tmp1) {
	    isplit[*nsplit] = j - 1;
	    ++(nsplit);
	    work[j - 1] = Zero;
	} else {
	    work[j - 1] = tmp1;
	    pivmin = max(pivmin, tmp1);
	}
    }
    isplit[*nsplit] = n;
    pivmin = pivmin * safemn;
//Compute interval and ATOLI
    if (irange == 3) {
//RANGE='I': Compute the INTEGERerval containing eigenvalues
//           IL through IU.
//Compute Gershgorin INTEGERerval for entire (split) matrix
//and use it as the initial INTEGERerval
	gu = d[1];
	gl = d[1];
	tmp1 = Zero;
	for (j = 0; j < n - 1; j++) {
	    tmp2 = sqrt(work[j]);
	    mtemp1 = gu, mtemp2 = d[j] + tmp1 + tmp2;
	    gu = max(mtemp1, mtemp2);
	    mtemp1 = gl, mtemp2 = d[j] - tmp1 - tmp2;
	    gl = min(mtemp1, mtemp2);
	    tmp1 = tmp2;
	}
	mtemp1 = gu, mtemp2 = d[n] + tmp1;
	gu = max(mtemp1, mtemp2);
	mtemp1 = gl, mtemp2 = d[n] - tmp1;
	gl = min(mtemp1, mtemp2);
	mtemp1 = abs(gl), mtemp2 = abs(gu);
	tnorm = max(mtemp1, mtemp2);
	gl = gl - tnorm * Fudge * ulp * n - pivmin * Fudge;
	gu = gu + tnorm * Fudge * ulp * n + pivmin * Fudge;
//Compute Iteration parameters
	itmax = (INTEGER) cast2double((log(tnorm + pivmin) - log(pivmin)) / log(Two) + 2.0);
	if (abstol <= Zero) {
	    atoli = ulp * tnorm;
	} else {
	    atoli = abstol;
	}
	work[n + 1] = gl;
	work[n + 2] = gl;
	work[n + 3] = gu;
	work[n + 4] = gu;
	work[n + 5] = gl;
	work[n + 6] = gu;
	iwork[1] = -1;
	iwork[2] = -1;
	iwork[3] = n + 1;
	iwork[4] = n + 1;
	iwork[5] = il - 1;
	iwork[6] = iu;
	Rlaebz(3, itmax, n, 2, 2, nb, atoli, rtoli, pivmin, &d[0], &e[0], &work[0], &iwork[5], &work[n + 1], &work[n + 5], &iout, &iwork[1], &w[1], &iblock[1], &iinfo);

	if (iwork[6] == iu) {
	    wl = work[n + 1];
	    wlu = work[n + 3];
	    nwl = iwork[1];
	    wu = work[n + 4];
	    wul = work[n + 2];
	    nwu = iwork[4];
	} else {
	    wl = work[n + 2];
	    wlu = work[n + 4];
	    nwl = iwork[2];
	    wu = work[n + 3];
	    wul = work[n + 1];
	    nwu = iwork[3];
	}
	if (nwl < 0 || nwl >= n || nwu < 1 || nwu > n) {
	    *info = 4;
	    return;
	}
    } else {
//RANGE='A' or 'V' -- Set ATOLI
	mtemp1 = abs(d[1]) + abs(e[1]), mtemp2 = abs(d[n]) + abs(e[n - 1]);
	tnorm = max(mtemp1, mtemp2);
	for (j = 2; j <= n - 1; j++) {
	    mtemp1 = tnorm, mtemp2 = abs(d[j]) + abs(e[j - 1]) + abs(e[j]);
	    tnorm = max(mtemp1, mtemp2);
	}
	if (abstol <= Zero) {
	    atoli = ulp * tnorm;
	} else {
	    atoli = abstol;
	}

	if (irange == 2) {
	    wl = vl;
	    wu = vu;
	} else {
	    wl = Zero;
	    wu = Zero;
	}
    }
//Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
//NWL accumulates the number of eigenvalues .le. WL,
//NWU accumulates the number of eigenvalues .le. WU
    m = 0;
    iend = 0;
    *info = 0;
    nwl = 0;
    nwu = 0;
    for (jb = 1; jb <= *nsplit; jb++) {
	ioff = iend;
	ibegin = ioff + 1;
	iend = isplit[jb];
	in = iend - ioff;
	if (in == 1) {
//Special Case -- IN=1
	    if (irange == 1 || wl >= d[ibegin] - pivmin) {
		++nwl;
	    }
	    if (irange == 1 || wu >= d[ibegin] - pivmin) {
		++nwu;
	    }
	    if (irange == 1 || (wl < d[ibegin] - pivmin && wu >= d[ibegin] - pivmin)) {
		++(*m);
		w[(*m)] = d[ibegin];
		iblock[(*m)] = jb;
	    }
	} else {
//General Case -- IN > 1
//Compute Gershgorin Interval
//and use it as the initial interval
	    gu = d[ibegin];
	    gl = d[ibegin];
	    tmp1 = Zero;
	    for (j = ibegin; j <= iend - 1; j++) {
		tmp2 = abs(e[j]);
		mtemp1 = gu, mtemp2 = d[j] + tmp1 + tmp2;
		gu = max(mtemp1, mtemp2);
		mtemp1 = gl, mtemp2 = d[j] - tmp1 - tmp2;
		gl = min(mtemp1, mtemp2);
		tmp1 = tmp2;
	    }
	    mtemp1 = gu, mtemp2 = d[iend] + tmp1;
	    gu = max(mtemp1, mtemp2);
	    mtemp1 = gl, mtemp2 = d[iend] - tmp1;
	    gl = min(mtemp1, mtemp2);
	    mtemp1 = abs(gl), mtemp2 = abs(gu);
	    bnorm = max(mtemp1, mtemp2);
	    gl = gl - bnorm * Fudge * ulp * in - pivmin * Fudge;
	    gu = gu + bnorm * Fudge * ulp * in + pivmin * Fudge;
//Compute ATOLI for the current submatrix
	    if (abstol <= Zero) {
		mtemp1 = abs(gl), mtemp2 = abs(gu);
		atoli = ulp * max(mtemp1, mtemp2);
	    } else {
		atoli = abstol;
	    }

	    if (irange > 1) {
		if (gu < wl) {
		    nwl = nwl + in;
		    nwu = nwu + in;
		    goto L70;
		}
		gl = max(gl, wl);
		gu = min(gu, wu);
		if (gl >= gu) {
		    goto L70;
		}
	    }
//Set Up Initial Mpackinterval
	    work[n + 1] = gl;
	    work[n + in + 1] = gu;
	    Rlaebz(1, 0, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin], &e[ibegin], &work[ibegin], &idumma, &work[n + 1],
		   &work[n + (in * 2) + 1], &im, &iwork[1], &w[(*m + 1)], &iblock[(*m + 1)], &iinfo);
	    nwl = nwl + iwork[1];
	    nwu = nwu + iwork[in + 1];
	    iwoff = (*m) - iwork[1];
//Compute Eigenvalues
	    itmax = (INTEGER) cast2double((log(gu - gl + pivmin) - log(pivmin)) / log(Two) + 2.0);
	    Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin], &e[ibegin], &work[ibegin], &idumma, &work[n + 1],
		   &work[n + (in * 2) + 1], &iout, &iwork[1], &w[(*m + 1)], &iblock[(*m + 1)], &iinfo);
//Copy Eigenvalues Mpackinto W and IBLOCK
//Use -JB for block number for unconverged eigenvalues.
	    for (j = 0; j < iout; j++) {
		tmp1 = (work[j + n] + work[j + in + n]) * Half;
//Flag non-convergence.
		if (j > iout - iinfo) {
		    ncnvrg = MTRUE;
		    ib = -jb;
		} else {
		    ib = jb;
		}
		for (je = iwork[j] + 1 + iwoff; je <= iwork[j + in] + iwoff; je++) {
		    w[je] = tmp1;
		    iblock[je] = ib;
		}
	    }
	    m = m + im;
	}
      L70:
	;
    }
//If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
//If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
    if (irange == 3) {
	im = 0;
	idiscl = il - 1 - nwl;
	idiscu = nwu - iu;
	if (idiscl > 0 || idiscu > 0) {
	    for (je = 1; je <= (*m); je++) {
		if (w[je] <= wlu && idiscl > 0) {
		    idiscl--;
		} else if (w[je] >= wul && idiscu > 0) {
		    idiscu--;
		} else {
		    im++;
		    w[im] = w[je];
		    iblock[im] = iblock[je];
		}
	    }
	    (*m) = im;
	}
	if (idiscl > 0 || idiscu > 0) {
//Code to deal with effects of bad arithmetic:
//Some low eigenvalues to be discarded are not in (WL,WLU],
//or high eigenvalues to be discarded are not in (WUL,WU]
//so just kill off the smallest IDISCL/largest IDISCU
//eigenvalues, by simply finding the smallest/largest
//eigenvalue(s).
//(If N(w) is monotone non-decreasing, this should never
//    happen.)
	    if (idiscl > 0) {
		wkill = wu;
		for (jdisc = 1; jdisc <= idiscl; jdisc++) {
		    iw = 0;
		    for (je = 1; je <= (*m); je++) {
			if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}
		    }
		    iblock[iw] = 0;
		}
	    }
	    if (idiscu > 0) {
		wkill = wl;
		for (jdisc = 1; jdisc <= idiscu; jdisc++) {
		    iw = 0;
		    for (je = 1; je <= (*m); je++) {
			if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}

		    }
		    iblock[iw] = 0;
		}
	    }
	    im = 0;
	    for (je = 1; je <= (*m); je++) {
		if (iblock[je] != 0) {
		    im++;
		    w[im] = w[je];
		    iblock[im] = iblock[je];
		}
	    }
	    (*m) = im;
	}
	if (idiscl < 0 || idiscu < 0) {
	    toofew = MTRUE;
	}
    }
//If ORDER='B', do nothing -- the eigenvalues are already sorted
//   by block.
//If ORDER='E', sort the eigenvalues from smallest to largest
    if (iorder == 1 && *nsplit > 1) {
	for (je = 1; je <= (*m - 1); je++) {
	    ie = 0;
	    tmp1 = w[je];
	    for (j = je + 1; j <= (*m); j++) {
		if (w[j] < tmp1) {
		    ie = j;
		    tmp1 = w[j];
		}
	    }
	    if (ie != 0) {
		itmp1 = iblock[ie];
		w[ie] = w[je];
		iblock[ie] = iblock[je];
		w[je] = tmp1;
		iblock[je] = itmp1;
	    }

	}
    }
    *info = 0;
    if (ncnvrg) {
	++(*info);
    }
    if (toofew) {
	*info += 2;
    }
    return;
}
