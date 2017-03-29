/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarrd.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlarrd(const char *range, const char *order, INTEGER n, REAL
	    vl, REAL vu, INTEGER il, INTEGER iu, REAL * gers,
	    REAL reltol, REAL * d, REAL * e, REAL * e2,
	    REAL pivmin, INTEGER nsplit, INTEGER * isplit, INTEGER * m,
	    REAL * w, REAL * werr, REAL * wl, REAL * wu, INTEGER * iblock, INTEGER * indexw, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, ib, ie, je, nb;
    REAL gl;
    INTEGER im, in;
    REAL gu;
    INTEGER iw, jee;
    REAL eps;
    INTEGER nwl;
    REAL wlu = 0.0, wul = 0.0;
    INTEGER nwu;
    REAL tmp1, tmp2;
    INTEGER iend, jblk, ioff, iout, itmp1, itmp2, jdisc;
    INTEGER iinfo;
    REAL atoli;
    INTEGER iwoff, itmax;
    REAL wkill, rtoli, uflow, tnorm;
    INTEGER ibegin;
    INTEGER irange, idiscl, idumma[1];
    REAL spdiam;
    INTEGER idiscu;
    INTEGER ncnvrg, toofew;
    REAL Zero = 0.0, Two = 2.0, Four = 4.0;
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
//Check for Errors
    if (irange <= 0) {
	*info = -1;
    } else if (!(Mlsame(order, "B") || Mlsame(order, "E"))) {
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
//Simplification:
    if (irange == 3 && il == 1 && iu == n) {
	irange = 1;
    }
//Get machine constants
    eps = Rlamch("P");
    uflow = Rlamch("U");
//Special Case when N=1
//Treat case of 1x1 matrix for quick return
    if (n == 1) {
	if (irange == 1 || (irange == 2 && d[1] > vl && d[1] <= vu) || (irange == 3 && il == 1 && iu == 1)) {
	    (*m) = 1;
	    w[1] = d[1];
//The computation error of the eigenvalue is zero
	    werr[1] = Zero;
	    iblock[1] = 1;
	    indexw[1] = 1;
	}
	return;
    }
//NB is the minimum vector length for vector bisection, or 0
//if only scalar is to be done.
    nb = iMlaenv(1, "Rstebz", " ", n, -1, -1, -1);
    if (nb <= 1) {
	nb = 0;
    }
//Find global spectral radius
    gl = d[1];
    gu = d[1];
    for (i = 0; i < n; i++) {
	mtemp1 = gl, mtemp2 = gers[(i * 2) - 1];
	gl = min(mtemp1, mtemp2);
	mtemp1 = gu, mtemp2 = gers[i * 2];
	gu = max(mtemp1, mtemp2);
    }
//Compute global Gerschgorin bounds and spectral diameter
    mtemp1 = abs(gl), mtemp2 = abs(gu);
    tnorm = max(mtemp1, mtemp2);
    gl = gl - tnorm * Two * eps * n - pivmin * Four;
    gu = gu + tnorm * Two * eps * n + pivmin * Four;
    spdiam = gu - gl;
//Input arguments for DLAEBZ:
//The relative tolerance.  An interval (a,b] lies within
//"relative tolerance" if  b-a < RELTOLmax(|a|,|b|),
    rtoli = reltol;
//Set the absolute tolerance for interval convergence to zero to force
//interval convergence based on relative size of the interval.
//This is dangerous because intervals might not converge when RELTOL is
//small. But at least a very small number should be selected so that for
//strongly graded matrices, the code can get relatively accurate
//eigenvalues.
    atoli = uflow * Four + pivmin * Four;
    if (irange == 3) {
//RANGE='I': Compute an interval containing eigenvalues
//IL through IU. The initial interval [GL,GU] from the global
//Gerschgorin bounds GL and GU is refined by DLAEBZ.
	itmax = (INTEGER) cast2double((log(tnorm + pivmin) - log(pivmin)) / log(Two)) + 2;
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
	Rlaebz(3, itmax, n, 2, 2, nb, atoli, rtoli, pivmin, &d[0], &e[0], &e2[1], &iwork[5], &work[n + 1], &work[n + 5]
	       , &iout, &iwork[1], &w[1], &iblock[1], &iinfo);
	if (iinfo != 0) {
	    *info = iinfo;
	    return;
	}
//On exit, output intervals may not be ordered by ascending negcount
	if (iwork[6] == iu) {
	    *wl = work[n + 1];
	    wlu = work[n + 3];
	    nwl = iwork[1];
	    *wu = work[n + 4];
	    wul = work[n + 2];
	    nwu = iwork[4];
	} else {
	    *wl = work[n + 2];
	    wlu = work[n + 4];
	    nwl = iwork[2];
	    *wu = work[n + 3];
	    wul = work[n + 1];
	    nwu = iwork[3];
	}
//On exit, the interval [WL, WLU] contains a value with negcount NWL,
//and [WUL, WU] contains a value with negcount NWU.
	if (nwl < 0 || nwl >= n || nwu < 1 || nwu > n) {
	    *info = 4;
	    return;
	}
    } else if (irange == 2) {
	*wl = vl;
	*wu = vu;
    } else if (irange == 1) {
	*wl = gl;
	*wu = gu;
    }
//Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
//NWL accumulates the number of eigenvalues .le. WL,
//NWU accumulates the number of eigenvalues .le. WU
    m = 0;
    iend = 0;
    *info = 0;
    nwl = 0;
    nwu = 0;

    for (jblk = 0; jblk <= nsplit; jblk++) {
	ioff = iend;
	ibegin = ioff + 1;
	iend = isplit[jblk];
	in = iend - ioff;
	if (in == 1) {
//1x1 block
	    if (*wl >= d[ibegin] - pivmin) {
		++nwl;
	    }
	    if (*wu >= d[ibegin] - pivmin) {
		++nwu;
	    }
	    if (irange == 1 || (*wl < d[ibegin] - pivmin && *wu >= d[ibegin] - pivmin)) {
		++(*m);
		w[(*m)] = d[ibegin];
		werr[(*m)] = Zero;
//The gap for a single block doesn't matter for the later
//algorithm and is assigned an arbitrary large value
		iblock[(*m)] = jblk;
		indexw[(*m)] = 1;
	    }
/*        Disabled 2x2 case because of a failure on the following matrix */
/*        RANGE = 'I', IL = IU = 4 */
/*          Original Tridiagonal, d = [ */
/*           -Zero150102010615740E+00 */
/*           -Zero849897989384260E+00 */
/*           -Zero128208148052635E-15 */
/*            Zero128257718286320E-15 */
/*          ]; */
/*          e = [ */
/*           -Zero357171383266986E+00 */
/*           -Zero180411241501588E-15 */
/*           -Zero175152352710251E-15 */
/*          ]; */
	} else {
/*           General Case - block of size IN >= 2 */
/*           Compute local Gerschgorin interval and use it as the initial */
/*           interval for DLAEBZ */
	    gu = d[ibegin];
	    gl = d[ibegin];
	    tmp1 = Zero;
	    for (j = ibegin; j <= iend; j++) {
		mtemp1 = gl, mtemp2 = gers[(j * 2) - 1];
		gl = min(mtemp1, mtemp2);
		mtemp1 = gu, mtemp2 = gers[j * 2];
		gu = max(mtemp1, mtemp2);
	    }
	    spdiam = gu - gl;
	    gl = gl - spdiam * Two * eps * in - pivmin * Two;
	    gu = gu + spdiam * Two * eps * in + pivmin * Two;
	    if (irange > 1) {
		if (gu < *wl) {
//the local block contains none of the wanted eigenvalues
		    nwl = nwl + in;
		    nwu = nwu + in;
		    goto L70;
		}
//refine search interval if possible, only range (WL,WU] matters
		gl = max(gl, *wl);
		gu = min(gu, *wu);
		if (gl >= gu) {
		    goto L70;
		}
	    }
//Find negcount of initial interval boundaries GL and GU
	    work[n + 1] = gl;
	    work[n + in + 1] = gu;
	    Rlaebz(1, 0, in, in, 1, nb, atoli, rtoli,
		   pivmin, &d[ibegin], &e[ibegin], &e2[ibegin], idumma, &work[n + 1], &work[n + (in << 1) + 1], &im, &iwork[1], &w[(*m) + 1], &iblock[(*m) + 1], &iinfo);
	    if (iinfo != 0) {
		*info = iinfo;
		return;
	    }
	    nwl = nwl + iwork[1];
	    nwu = nwu + iwork[in + 1];
	    iwoff = (*m) - iwork[1];
//Compute Eigenvalues
	    itmax = (INTEGER) cast2double((log(gu - gl + pivmin) - log(pivmin)) / log(Two)) + 2;
	    Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli,
		   pivmin, &d[ibegin], &e[ibegin], &e2[ibegin], idumma, &work[n + 1], &work[n + (in << 1) + 1], &iout, &iwork[1], &w[(*m) + 1], &iblock[(*m) + 1], &iinfo);
	    if (iinfo != 0) {
		*info = iinfo;
		return;
	    }
//Copy eigenvalues into W and IBLOCK
//Use -JBLK for block number for unconverged eigenvalues.
//Loop over the number of output intervals from DLAEBZ
	    for (j = 0; j < iout; j++) {
//eigenvalue approximation is middle point of interval
		tmp1 = (work[j + n] + work[j + in + n]) * .5;
//semi length of error interval
		tmp2 = abs(work[j + n] - work[j + in + n]) * .5;
		if (j > iout - iinfo) {
//Flag non-convergence.
		    ncnvrg = MTRUE;
		    ib = -jblk;
		} else {
		    ib = jblk;
		}
		for (je = iwork[j] + 1 + iwoff; je <= iwork[j + in] + iwoff; je++) {
		    w[je] = tmp1;
		    werr[je] = tmp2;
		    indexw[je] = je - iwoff;
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
	idiscl = il - 1 - nwl;
	idiscu = nwu - iu;
	if (idiscl > 0) {
	    im = 0;
	    for (je = 1; je <= (*m); je++) {
//Remove some of the smallest eigenvalues from the left so that
//at the end IDISCL =Zero Move all eigenvalues up to the left.
		if (w[je] <= wlu && idiscl > 0) {
		    idiscl--;
		} else {
		    im++;
		    w[im] = w[je];
		    werr[im] = werr[je];
		    indexw[im] = indexw[je];
		    iblock[im] = iblock[je];
		}
	    }
	    (*m) = im;
	}
	if (idiscu > 0) {
//Remove some of the largest eigenvalues from the right so that
//at the end IDISCU =Zero Move all eigenvalues up to the left.
	    im = (*m) + 1;
	    for (je = (*m); je >= 1; je--) {
		if (w[je] >= wul && idiscu > 0) {
		    idiscu--;
		} else {
		    im--;
		    w[im] = w[je];
		    werr[im] = werr[je];
		    indexw[im] = indexw[je];
		    iblock[im] = iblock[je];
		}
	    }
	    jee = 0;
	    for (je = im; je <= (*m); je++) {
		jee++;
		w[jee] = w[je];
		werr[jee] = werr[je];
		indexw[jee] = indexw[je];
		iblock[jee] = iblock[je];
	    }
	    m = m - im + 1;
	}
	if (idiscl > 0 || idiscu > 0) {
//Code to deal with effects of bad arithmetic. (If N(w) is
//monotone non-decreasing, this should never happen.)
//Some low eigenvalues to be discarded are not in (WL,WLU],
//or high eigenvalues to be discarded are not in (WUL,WU]
//so just kill off the smallest IDISCL/largest IDISCU
//eigenvalues, by marking the corresponding IBLOCK = 0
	    if (idiscl > 0) {
		wkill = *wu;
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
		wkill = *wl;
		for (jdisc = 1; jdisc <= idiscu; jdisc++) {
		    iw = 0;
		    for (je = 1; je <= (*m); je++) {
			if (iblock[je] != 0 && (w[je] >= wkill || iw == 0)) {
			    iw = je;
			    wkill = w[je];
			}

		    }
		    iblock[iw] = 0;
		}
	    }
//Now erase all eigenvalues with IBLOCK set to zero
	    im = 0;
	    for (je = 1; je <= (*m); je++) {
		if (iblock[je] != 0) {
		    im++;
		    w[im] = w[je];
		    werr[im] = werr[je];
		    indexw[im] = indexw[je];
		    iblock[im] = iblock[je];
		}
	    }
	    (*m) = im;
	}
	if (idiscl < 0 || idiscu < 0) {
	    toofew = MTRUE;
	}
    }
    if ((irange == 1 && (*m) != n) || (irange == 3 && (*m) != iu - il + 1)) {
	toofew = MTRUE;
    }
//If ORDER='B', do nothing the eigenvalues are already sorted by
//   block.
//If ORDER='E', sort the eigenvalues from smallest to largest
    if (Mlsame(order, "E") && nsplit > 1) {
	for (je = 1; je <= (*m) - 1; je++) {
	    ie = 0;
	    tmp1 = w[je];
	    for (j = je + 1; j <= (*m); j++) {
		if (w[j] < tmp1) {
		    ie = j;
		    tmp1 = w[j];
		}
	    }
	    if (ie != 0) {
		tmp2 = werr[ie];
		itmp1 = iblock[ie];
		itmp2 = indexw[ie];
		w[ie] = w[je];
		werr[ie] = werr[je];
		iblock[ie] = iblock[je];
		indexw[ie] = indexw[je];
		w[je] = tmp1;
		werr[je] = tmp2;
		iblock[je] = itmp1;
		indexw[je] = itmp2;
	    }
	}
    }
    *info = 0;
    if (ncnvrg) {
	++(*info);
    }
    if (toofew) {
	*info = *info + 2;
    }
    return;
}
