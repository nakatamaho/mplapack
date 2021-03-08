/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarre.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlarre(const char *range, INTEGER n, REAL * vl, REAL * vu, INTEGER il, INTEGER iu, REAL * d, REAL
	    * e, REAL * e2, REAL rtol1, REAL rtol2, REAL spltol, INTEGER * nsplit, INTEGER * isplit,
	    INTEGER * m, REAL * w, REAL * werr, REAL * wgap, INTEGER * iblock, INTEGER * indexw, REAL * gers, REAL * pivmin, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j;
    REAL s1, s2;
    INTEGER mb = 0;
    REAL gl;
    INTEGER in, mm;
    REAL gu;
    INTEGER cnt;
    REAL eps, tau, tmp, rtl;
    INTEGER cnt1, cnt2;
    REAL tmp1, eabs;
    INTEGER iend, jblk;
    REAL eold;
    INTEGER indl = 0;
    REAL dmax, emax;
    INTEGER wend = 0, idum, indu = 0;
    REAL rtol;
    INTEGER iseed[4];
    REAL avgap, sigma;
    INTEGER iinfo;
    INTEGER norep;
    INTEGER ibegin;
    INTEGER forceb;
    INTEGER irange = 0;
    REAL sgndef;
    INTEGER wbegin;
    REAL safmin, spdiam;
    INTEGER usedqd = 0;
    REAL clwdth, isleft;
    REAL isrght, bsrtol, dpivot;
    REAL Zero = 0.0, Half = .5, One = 1.0, Two = 2.0, Four = 4.0, Eight = 8.0, Fourth = .25, Hundrd = 100.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    if (Mlsame(range, "A")) {
	irange = 1;
    } else if (Mlsame(range, "V")) {
	irange = 3;
    } else if (Mlsame(range, "I")) {
	irange = 2;
    }
    *m = 0;
//Get machine constants
    safmin = Rlamch("S");
    eps = Rlamch("P");
//Set parameters
    rtl = sqrt(eps);
    bsrtol = sqrt(eps);
//Treat case of 1x1 matrix for quick return
    if (n == 1) {
	if (irange == 1 || (irange == 3 && d[1] > *vl && d[1] <= *vu) || (irange == 2 && il == 1 && iu == 1)) {
	    (*m) = 1;
	    w[1] = d[1];
//The computation error of the eigenvalue is zero
	    werr[1] = Zero;
	    wgap[1] = Zero;
	    iblock[1] = 1;
	    indexw[1] = 1;
	    gers[1] = d[1];
	    gers[2] = d[1];
	}
//store the shift for the initial RRR, which is zero in this case
	e[1] = Zero;
	return;
    }
//General case: tridiagonal matrix of order > 1
//Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
//Compute maximum off-diagonal entry and pivmin.
    gl = d[1];
    gu = d[1];
    eold = Zero;
    emax = Zero;
    e[n] = Zero;
    for (i = 0; i < n; i++) {
	werr[i] = Zero;
	wgap[i] = Zero;
	eabs = abs(e[i]);
	if (eabs >= emax) {
	    emax = eabs;
	}
	tmp1 = eabs + eold;
	gers[i * 2 - 1] = d[i] - tmp1;
	mtemp1 = gl, mtemp2 = gers[i * 2 - 1];
	gl = min(mtemp1, mtemp2);
	gers[i * 2] = d[i] + tmp1;
	mtemp1 = gu, mtemp2 = gers[i * 2];
	gu = max(mtemp1, mtemp2);
	eold = eabs;
    }
//The minimum pivot allowed in the Sturm sequence for T
    mtemp1 = One, mtemp2 = emax * emax;
    *pivmin = safmin * max(mtemp1, mtemp2);
//Compute spectral diameter. The Gerschgorin bounds give an
//estimate that is wrong by at most a factor of SQRT(2)
    spdiam = gu - gl;
//Compute splitting points
    Rlarra(n, &d[0], &e[0], &e2[1], spltol, spdiam, nsplit, &isplit[1], &iinfo);
//Can force use of bisection instead of faster DQDS.
//Option left in the code for future multisection work.
    forceb = MFALSE;
    if (irange == 1 && !forceb) {
//Set interval [VL,VU] that contains all eigenvalues
	*vl = gl;
	*vu = gu;
    } else {
//We call DLARRD to find crude approximations to the eigenvalues
//in the desired range. In case IRANGE = INDRNG, we also obtain the
//interval (VL,VU] that contains all the wanted eigenvalues.
//An interval [LEFT,RIGHT] has converged if
//RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT))
//DLARRD needs a WORK of size 4*N, IWORK of size 3*N
	Rlarrd(range, "B", n, (*vl), (*vu), il, iu, &gers[1], bsrtol, &d[0], &e[0], &e2[1], *pivmin, *nsplit, &isplit[1], &mm,
	       &w[1], &werr[1], vl, vu, &iblock[1], &indexw[1], &work[0], &iwork[1], &iinfo);
	if (iinfo != 0) {
	    *info = -1;
	    return;
	}
//Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0
	for (i = mm + 1; i <= n; i++) {
	    w[i] = Zero;
	    werr[i] = Zero;
	    iblock[i] = 0;
	    indexw[i] = 0;
	}
    }
//Loop over unreduced blocks
    ibegin = 1;
    wbegin = 1;
    for (jblk = 0; jblk <= (*nsplit); jblk++) {
	iend = isplit[jblk];
	in = iend - ibegin + 1;
//1 X 1 block
	if (in == 1) {
	    if (irange == 1 || (irange == 3 && d[ibegin] > *vl && d[ibegin] <= *vu) || (irange == 2 && iblock[wbegin] == jblk)) {
		++(*m);
		w[(*m)] = d[ibegin];
		werr[(*m)] = Zero;
//The gap for a single block doesn't matter for the later
//algorithm and is assigned an arbitrary large value
		wgap[(*m)] = Zero;
		iblock[(*m)] = jblk;
		indexw[(*m)] = 1;
		++wbegin;
	    }
//E( IEND ) holds the shift for the initial RRR
	    e[iend] = Zero;
	    ibegin = iend + 1;
	    goto L170;
	}
//Blocks of size larger than 1x1
//E( IEND ) will hold the shift for the initial RRR, for now set it =0
	e[iend] = Zero;
//Find local outer bounds GL,GU for the block
	gl = d[ibegin];
	gu = d[ibegin];
	for (i = ibegin; i <= iend; i++) {
	    mtemp1 = gers[i * 2 - 1], mtemp2 = gl;
	    gl = min(mtemp1, mtemp2);
	    mtemp1 = gers[i * 2], mtemp2 = gu;
	    gu = max(mtemp1, mtemp2);
	}
	spdiam = gu - gl;
	if (!(irange == 1 && !forceb)) {
//Count the number of eigenvalues in the current block.
	    mb = 0;
	    for (i = wbegin; i <= mm; i++) {
		if (iblock[i] == jblk) {
		    mb++;
		} else {
		    goto L21;
		}
	    }
	  L21:
	    if (mb == 0) {
//No eigenvalue in the current block lies in the desired range
//E( IEND ) holds the shift for the initial RRR
		e[iend] = Zero;
		ibegin = iend + 1;
		goto L170;
	    } else {
//Decide whether dqds or bisection is more efficient
		usedqd = mb > in * Half && !forceb;
		wend = wbegin + mb - 1;
//Calculate gaps for the current block
//In later stages, when representations for individual
//eigenvalues are different, we use SIGMA = E( IEND ).
		sigma = Zero;
		for (i = wbegin; i <= wend - 1; i++) {
		    mtemp1 = Zero, mtemp2 = w[i + 1] - werr[i + 1] - (w[i] + werr[i]);
		    wgap[i] = max(mtemp1, mtemp2);
		}
		mtemp1 = Zero, mtemp2 = *vu - sigma - (w[wend] + werr[wend]);
		wgap[wend] = max(mtemp1, mtemp2);
//Find local index of the first and last desired evalue.
		indl = indexw[wbegin];
		indu = indexw[wend];
	    }
	}
	if ((irange == 1 && !forceb) || usedqd) {
//Case of DQDS
//Find approximations to the extremal eigenvalues of the block
	    Rlarrk(in, 1, gl, gu, &d[ibegin], &e2[ibegin], (*pivmin), rtl, &tmp, &tmp1, &iinfo);
	    if (iinfo != 0) {
		*info = -1;
		return;
	    }
	    mtemp1 = gl, mtemp2 = tmp - tmp1 - eps * 100 * abs(tmp - tmp1);
	    isleft = max(mtemp1, mtemp2);
	    Rlarrk(in, in, gl, gu, &d[ibegin], &e2[ibegin], (*pivmin), rtl, &tmp, &tmp1, &iinfo);
	    if (iinfo != 0) {
		*info = -1;
		return;
	    }
	    mtemp1 = gu, mtemp2 = tmp + tmp1 + eps * Hundrd * abs(tmp + tmp1);
	    isrght = min(mtemp1, mtemp2);
//Improve the estimate of the spectral diameter
	    spdiam = isrght - isleft;
	} else {
//Case of bisection
//Find approximations to the wanted extremal eigenvalues
	    mtemp1 = gl, mtemp2 = w[wbegin] - werr[wbegin] - eps * Hundrd * abs(w[wbegin] - werr[wbegin]);
	    isleft = max(mtemp1, mtemp2);
	    mtemp1 = gu, mtemp2 = w[wend] + werr[wend] + eps * Hundrd * abs(w[wend] + werr[wend]);
	    isrght = min(mtemp1, mtemp2);
	}
//Decide whether the base representation for the current block
//L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
//should be on the left or the right end of the current block.
//The strategy is to shift to the end which is "more populated"
//Furthermore, decide whether to use DQDS for the computation of
//the eigenvalue approximations at the end of DLARRE or bisection.
//dqds is chosen if all eigenvalues are desired or the number of
//eigenvalues to be computed is large compared to the blocksize.
	if (irange == 1 && !forceb) {
//If all the eigenvalues have to be computed, we use dqd
	    usedqd = MTRUE;
//INDL is the local index of the first eigenvalue to compute
	    indl = 0;
	    indu = in;
//MB =  number of eigenvalues to compute
	    mb = in;
	    wend = wbegin + mb - 1;
//Define 1/4 and 3/4 points of the spectrum
	    s1 = isleft + spdiam * Fourth;
	    s2 = isrght - spdiam * Fourth;
	} else {
//DLARRD has computed IBLOCK and INDEXW for each eigenvalue
//approximation.
//choose sigma
	    if (usedqd) {
		s1 = isleft + spdiam * Fourth;
		s2 = isrght - spdiam * Fourth;
	    } else {
		tmp = min(isrght, *vu) - max(isleft, *vl);
		s1 = max(isleft, *vl) + tmp * Fourth;
		s2 = min(isrght, *vu) - tmp * Fourth;
	    }
	}
//Compute the negcount at the 1/4 and 3/4 points
	if (mb > 1) {
	    Rlarrc("T", in, s1, s2, &d[ibegin], &e[ibegin], *pivmin, &cnt, &cnt1, &cnt2, &iinfo);
	}
	if (mb == 1) {
	    sigma = gl;
	    sgndef = One;
	} else if (cnt1 - indl >= indu - cnt2) {
	    if (irange == 1 && !forceb) {
		sigma = max(isleft, gl);
	    } else if (usedqd) {
//use Gerschgorin bound as shift to get pos def matrix
//for dqds
		sigma = isleft;
	    } else {
//use approximation of the first desired eigenvalue of the
//block as shift
		sigma = max(isleft, *vl);
	    }
	    sgndef = One;
	} else {
	    if (irange == 1 && !forceb) {
		sigma = min(isrght, gu);
	    } else if (usedqd) {
//use Gerschgorin bound as shift to get neg def matrix
//for dqds
		sigma = isrght;
	    } else {
//use approximation of the first desired eigenvalue of the
//block as shift
		sigma = min(isrght, *vu);
	    }
	    sgndef = -One;
	}
//An initial SIGMA has been chosen that will be used for computing
//T - SIGMA I = L D L^T
//Define the increment TAU of the shift in case the initial shift
//needs to be refined to obtain a factorization with not too much
//element growth.
	if (usedqd) {
//The initial SIGMA was to the outer end of the spectrum
//the matrix is definite and we need not retreat.
	    tau = spdiam * eps * n + (*pivmin) * Two;
	} else {
	    if (mb > 1) {
		clwdth = w[wend] + werr[wend] - w[wbegin] - werr[wbegin];
		avgap = clwdth / abs(wend - wbegin);
		if (sgndef == One) {
		    mtemp1 = wgap[wbegin];
		    tau = max(mtemp1, avgap) * Half;
		    mtemp1 = tau, mtemp2 = werr[wbegin];
		    tau = max(mtemp1, mtemp2);
		} else {
		    mtemp1 = wgap[wend - 1];
		    tau = max(mtemp1, avgap) * Half;
		    mtemp1 = tau, mtemp2 = werr[wend];
		    tau = max(mtemp1, mtemp2);
		}
	    } else {
		tau = werr[wbegin];
	    }
	}
	for (idum = 1; idum <= 6; idum++) {
//Compute L D L^T factorization of tridiagonal matrix T - sigma I.
//Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of
//pivots in WORK(2*IN+1:3*IN)
	    dpivot = d[ibegin] - sigma;
	    work[1] = dpivot;
	    dmax = abs(work[1]);
	    j = ibegin;
	    for (i = 0; i < in - 1; i++) {
		work[(in << 1) + i] = One / work[i];
		tmp = e[j] * work[(in << 1) + i];
		work[in + i] = tmp;
		dpivot = d[j + 1] - sigma - tmp * e[j];
		work[i + 1] = dpivot;
		mtemp1 = dmax, mtemp2 = abs(dpivot);
		dmax = max(mtemp1, mtemp2);
		j++;
	    }
//check for element growth
	    if (dmax > spdiam * 64) {
		norep = MTRUE;
	    } else {
		norep = MFALSE;
	    }
	    if (usedqd && !norep) {
//Ensure the definiteness of the representation
//All entries of D (of L D L^T) must have the same sign
		for (i = 0; i < in; i++) {
		    tmp = sgndef * work[i];
		    if (tmp < Zero) {
			norep = MTRUE;
		    }
		}
	    }
	    if (norep) {
//Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
//shift which makes the matrix definite. So we should end up
//here really only in the case of IRANGE = VALRNG or INDRNG.
		if (idum == 5) {
		    if (sgndef == One) {
//The fudged Gerschgorin shift should succeed
			sigma = gl - spdiam * Two * eps * n - (*pivmin) * Four;
		    } else {
			sigma = gu + spdiam * Two * eps * n + (*pivmin) * Four;
		    }
		} else {
		    sigma = sigma - sgndef * tau;
		    tau = tau * Two;
		}
	    } else {
//an initial RRR is found
		goto L83;
	    }
	}
//if the program reaches this point, no base representation could be
//found in MAXTRY iterations.
	*info = 2;
	return;
      L83:
//At this point, we have found an initial base representation
//T - SIGMA I = L D L^T with not too much element growth.
//Store the shift.
	e[iend] = sigma;
//Store D and L.
	Rcopy(in, &work[0], 1, &d[ibegin], 1);
	Rcopy(in - 1, &work[in + 1], 1, &e[ibegin], 1);
	if (mb > 1) {
//Perturb each entry of the base representation by a small
//(but random) relative amount to overcome difficulties with
//glued matrices.
	    for (i = 1; i < 4; i++) {
		iseed[i - 1] = 1;
	    }
	    Rlarnv(2, iseed, (in * 2) - 1, &work[0]);
	    for (i = 0; i < in - 1; i++) {
		d[ibegin + i - 1] = d[ibegin + i - 1] * (eps * Eight * work[i] + One);
		e[ibegin + i - 1] = e[ibegin + i - 1] * (eps * Eight * work[in + i] + One);
	    }
	    d[iend] = d[iend] * (eps * Four * work[in] + One);
	}
//Don't update the Gerschgorin intervals because keeping track
//of the updates would be too much work in DLARRV.
//We update W instead and use it to locate the proper Gerschgorin
//intervals.
//Compute the required eigenvalues of L D L' by bisection or dqds
	if (!usedqd) {
//If DLARRD has been used, shift the eigenvalue approximations
//according to their representation. This is necessary for
//a uniform DLARRV since dqds computes eigenvalues of the
//shifted representation. In DLARRV, W will always hold the
//Unshifted eigenvalue approximation.
	    for (j = wbegin; j <= wend; j++) {
		w[j] = w[j] - sigma;
		werr[j] = werr[j] + abs(w[j]) * eps;
	    }
//call DLARRB to reduce eigenvalue error of the approximations
//from DLARRD
	    for (i = ibegin; i <= iend - 1; i++) {
		work[i] = d[i] * (e[i] * e[i]);
	    }
//use bisection to find EV from INDL to INDU
	    Rlarrb(in, &d[ibegin], &work[ibegin], indl, indu, rtol1,
		   rtol2, indl - 1, &w[wbegin], &wgap[wbegin], &werr[wbegin], &work[n * 2 + 1], &iwork[1], (*pivmin), spdiam, in, &iinfo);
	    if (iinfo != 0) {
		*info = -4;
		return;
	    }
//DLARRB computes all gaps correctly except for the last one
//Record distance to VU/GU
	    mtemp1 = Zero, mtemp2 = *vu - sigma - (w[wend] + werr[wend]);
	    wgap[wend] = max(mtemp1, mtemp2);
	    for (i = indl; i <= indu; i++) {
		++(*m);
		iblock[(*m)] = jblk;
		indexw[(*m)] = i;
	    }
	} else {
//Call dqds to get all eigs (and then possibly delete unwanted
//eigenvalues).
//Note that dqds finds the eigenvalues of the L D L^T representation
//of T to high relative accuracy. High relative accuracy
//might be lost when the shift of the RRR is subtracted to obtain
//the eigenvalues of T. However, T is not guaranteed to define its
//eigenvalues to high relative accuracy anyway.
//Set RTOL to the order of the tolerance used in DLASQ2
//This is an ESTIMATED error, the worst case bound is 4*N*EPS
//which is usually too large and requires unnecessary work to be
//done by bisection when computing the eigenvectors
	    rtol = log(in) * Four * eps;
	    j = ibegin;
	    for (i = 0; i < in - 1; i++) {
		work[i * 2 - 1] = abs(d[j]);
		work[i * 2] = e[j] * e[j] * work[i * 2 - 1];
		j++;

	    }
	    work[in * 2 - 1] = abs(d[iend]);
	    work[in * 2] = Zero;
	    Rlasq2(in, &work[0], &iinfo);
	    if (iinfo != 0) {
//If IINFO = -5 then an index is part of a tight cluster
//and should be changed. The index is in IWORK(1) and the
//gap is in WORK(N+1)
		*info = -5;
		return;
	    } else {
//Test that all eigenvalues are positive as expected
		for (i = 0; i < in; i++) {
		    if (work[i] < Zero) {
			*info = -6;
			return;
		    }
		}
	    }
	    if (sgndef > Zero) {
		for (i = indl; i <= indu; i++) {
		    ++(*m);
		    w[(*m)] = work[in - i + 1];
		    iblock[(*m)] = jblk;
		    indexw[(*m)] = i;
		}
	    } else {
		for (i = indl; i <= indu; i++) {
		    ++(*m);
		    w[(*m)] = -work[i];
		    iblock[(*m)] = jblk;
		    indexw[(*m)] = i;
		}
	    }
	    for (i = (*m) - mb + 1; i <= (*m); i++) {
//the value of RTOL below should be the tolerance in DLASQ2
		werr[i] = rtol * abs(w[i]);
	    }
	    for (i = (*m) - mb + 1; i <= (*m) - 1; i++) {
//compute the right gap between the intervals
		mtemp1 = Zero, mtemp2 = w[i + 1] - werr[i + 1] - (w[i] + werr[i]);
		wgap[i] = max(mtemp1, mtemp2);
	    }
	    mtemp1 = Zero, mtemp2 = *vu - sigma - (w[(*m)] + werr[(*m)]);
	    wgap[(*m)] = max(mtemp1, mtemp2);
	}
//proceed with next block
	ibegin = iend + 1;
	wbegin = wend + 1;
      L170:;
    }
    return;
}
