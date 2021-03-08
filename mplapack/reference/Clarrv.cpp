/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarrv.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clarrv(INTEGER n, REAL vl, REAL vu,
	    REAL * d, REAL * l, REAL pivmin, INTEGER * isplit,
	    INTEGER m, INTEGER dol, INTEGER dou, REAL minrgp,
	    REAL rtol1, REAL rtol2, REAL * w, REAL * werr,
	    REAL * wgap, INTEGER * iblock, INTEGER * indexw, REAL * gers, COMPLEX * z, INTEGER ldz, INTEGER * isuppz, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER minwsize, i, j, k, p, q, miniwsize, ii;
    REAL gl;
    INTEGER im, in;
    REAL gu, gap, eps, tau, tol, tmp;
    INTEGER zto;
    REAL ztz;
    INTEGER iend, jblk;
    REAL lgap;
    INTEGER done;
    REAL rgap, left;
    INTEGER wend, iter;
    REAL bstw = 0;
    INTEGER itmp1, indld;
    REAL fudge;
    INTEGER idone;
    REAL sigma;
    INTEGER iinfo, iindr;
    REAL resid;
    INTEGER eskip;
    REAL right;
    INTEGER nclus, zfrom;
    REAL rqtol;
    INTEGER iindc1, iindc2, indin1, indin2;
    INTEGER stp2ii;
    INTEGER ibegin, indeig;
    INTEGER needbs;
    INTEGER indlld;
    REAL sgndef, mingma;
    INTEGER oldien, oldncl, wbegin;
    REAL spdiam;
    INTEGER negcnt;
    INTEGER oldcls;
    REAL savgap = 0;
    INTEGER ndepth;
    REAL ssigma;
    LOGICAL usedbs;
    INTEGER iindwk, offset;
    REAL gaptol;
    INTEGER newcls, oldfst, indwrk, windex, oldlst;
    INTEGER usedrq;
    INTEGER newfst, newftt, parity, windmn, windpl, isupmn, newlst, zusedl;
    REAL bstres = 0;
    INTEGER newsiz, zusedu, zusedw;
    REAL nrminv, lambda;
    INTEGER tryrqc;
    INTEGER isupmx;
    REAL rqcorr;
    REAL Zero = 0.0, Half = 0.5, One = 1.0, Two = 2.0, Three = 3.0, Four = 4.0;
    REAL mtemp1, mtemp2;

//The first N entries of WORK are reserved for the eigenvalues
    indld = n + 1;
    indlld = n * 2 + 1;
    indin1 = n * 3 + 1;
    indin2 = n * 4 + 1;
    indwrk = n * 5 + 1;
    minwsize = n * 12;
    for (i = 0; i < minwsize; i++) {
	work[i] = Zero;
    }
//IWORK(IINDR+1:IINDR+N) hold the twist indices R for the
//factorization used to compute the FP vector
    iindr = 0;
//IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
//layer and the one above.
    iindc1 = n;
    iindc2 = n * 2;
    iindwk = n * 3 + 1;
    miniwsize = n * 7;
    for (i = 0; i < miniwsize; i++) {
	iwork[i] = 0;
    }
    zusedl = 0;
    if (dol > 1) {
//Set lower bound for use of Z
	zusedl = dol - 1;
    }
    zusedu = m;
    if (dou < m) {
//Set lower bound for use of Z
	zusedu = dou + 1;
    }
//The width of the part of Z that is used
    zusedw = zusedu - zusedl + 1;
    Claset("Full", n, zusedw, Zero, Zero, &z[zusedl * ldz + 1], ldz);
    eps = Rlamch("Precision");
    rqtol = eps * Two;
//Set expert flags for standard code.
    tryrqc = MTRUE;
    if (dol == 1 && dou == m) {
    } else {
//Only selected eigenpairs are computed. Since the other evalues
//are not refined by RQ iteration, bisection has to compute to full
//accuracy.
	rtol1 = eps * Four;
	rtol2 = eps * Four;
    }
//The entries WBEGIN:WEND in W, WERR, WGAP correspond to the
//desired eigenvalues. The support of the nonzero eigenvector
//entries is contained in the interval IBEGIN:IEND.
//Remark that if k eigenpairs are desired, then the eigenvectors
//are stored in k contiguous columns of Z.
//DONE is the number of eigenvectors already computed
    done = 0;
    ibegin = 1;
    wbegin = 1;
    for (jblk = 0; jblk <= iblock[m]; jblk++) {
	iend = isplit[jblk];
	sigma = l[iend];
//Find the eigenvectors of the submatrix indexed IBEGIN
//through IEND.
	wend = wbegin - 1;
      L15:
	if (wend < m) {
	    if (iblock[wend + 1] == jblk) {
		++wend;
		goto L15;
	    }
	}
	if (wend < wbegin) {
	    ibegin = iend + 1;
	    goto L170;
	} else if (wend < dol || wbegin > dou) {
	    ibegin = iend + 1;
	    wbegin = wend + 1;
	    goto L170;
	}
//Find local spectral diameter of the block
	gl = gers[ibegin * 2 - 1];
	gu = gers[ibegin * 2];
	for (i = ibegin + 1; i <= iend; i++) {
	    mtemp1 = gers[(i * 2) - 1], mtemp2 = gl;
	    gl = min(mtemp1, mtemp2);
	    mtemp1 = gers[i * 2], mtemp2 = gu;
	    gu = max(mtemp1, mtemp2);
	}
	spdiam = gu - gl;
//OLDIEN is the last index of the previous block
	oldien = ibegin - 1;
//Calculate the size of the current block
	in = iend - ibegin + 1;
//The number of eigenvalues in the current block
	im = wend - wbegin + 1;
//This is for a 1x1 block
	if (ibegin == iend) {
	    ++done;
	    z[ibegin + wbegin * ldz] = One;
	    isuppz[(wbegin * 2) - 1] = ibegin;
	    isuppz[wbegin * 2] = ibegin;
	    w[wbegin] = w[wbegin] + sigma;
	    work[wbegin] = w[wbegin];
	    ibegin = iend + 1;
	    ++wbegin;
	    goto L170;
	}
//The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND)
//Note that these can be approximations, in this case, the corresp.
//entries of WERR give the size of the uncertainty interval.
//The eigenvalue approximations will be refined when necessary as
//high relative accuracy is required for the computation of the
//corresponding eigenvectors.
	Rcopy(im, &w[wbegin], 1, &work[wbegin], 1);
//We store in W the eigenvalue approximations w.r.t. the original
//matrix T.
	for (i = 0; i < im; i++) {
	    w[wbegin + i - 1] = w[wbegin + i - 1] + sigma;
	}
//NDEPTH is the current depth of the representation tree
	ndepth = 0;
//PARITY is either 1 or 0
	parity = 1;
//NCLUS is the number of clusters for the next level of the
//representation tree, we start with NCLUS = 1 for the root
	nclus = 1;
	iwork[iindc1 + 1] = 1;
	iwork[iindc1 + 2] = im;
//IDONE is the number of eigenvectors already computed in the current
//block
	idone = 0;
//loop while( IDONE.LT.IM )
//generate the representation tree for the current block and
//compute the eigenvectors
      L40:
	if (idone < im) {
//This is a crude protection against infinitely deep trees
	    if (ndepth > m) {
		*info = -2;
		return;
	    }
//breadth first processing of the current level of the representation
//tree: OLDNCL = number of clusters on current level
	    oldncl = nclus;
//reset NCLUS to count the number of child clusters
	    nclus = 0;
	    parity = 1 - parity;
	    if (parity == 0) {
		oldcls = iindc1;
		newcls = iindc2;
	    } else {
		oldcls = iindc2;
		newcls = iindc1;
	    }
//Process the clusters on the current level
	    for (i = 0; i < oldncl; i++) {
		j = oldcls + (i * 2);
//OLDFST, OLDLST = first, last index of current cluster.
//                 cluster indices start with 1 and are relative
//                 to WBEGIN when accessing W, WGAP, WERR, Z
		oldfst = iwork[j - 1];
		oldlst = iwork[j];
		if (ndepth > 0) {
//Retrieve relatively robust representation (RRR) of cluster
//that has been computed at the previous level
//The RRR is stored in Z and overwritten once the eigenvectors
//have been computed or when the cluster is refined
		    if (dol == 1 && dou == m) {
//Get representation from location of the leftmost evalue
//of the cluster
			j = wbegin + oldfst - 1;
		    } else {
			if (wbegin + oldfst - 1 < dol) {
//Get representation from the left end of Z array
			    j = dol - 1;
			} else if (wbegin + oldfst - 1 > dou) {
//Get representation from the right end of Z array
			    j = dou;
			} else {
			    j = wbegin + oldfst - 1;
			}
		    }
		    for (k = 0; k < in - 1; k++) {
			d[ibegin + k - 1] = z[ibegin + k - 1 + j * ldz].real();
			l[ibegin + k - 1] = z[ibegin + k - 1 + (j + 1) * ldz].real();
		    }
		    d[iend] = z[iend + j * ldz].real();
		    sigma = z[iend + (j + 1) * ldz].real();
//Set the corresponding entries in Z to zero
		    Claset("Full", in, 2, Zero, Zero, &z[ibegin + j * ldz], ldz);
		}
//Compute DL and DLL of current RRR
		for (j = ibegin; j <= iend - 1; j++) {
		    tmp = d[j] * l[j];
		    work[indld - 1 + j] = tmp;
		    work[indlld - 1 + j] = tmp * l[j];

		}
		if (ndepth > 0) {
//P and Q are index of the first and last eigenvalue to compute
//within the current block
		    p = indexw[wbegin - 1 + oldfst];
		    q = indexw[wbegin - 1 + oldlst];
//Offset for the arrays WORK, WGAP and WERR, i.e., th P-OFFSET
//thru' Q-OFFSET elements of these arrays are to be used.
// OFFSET = P-OLDFST
		    offset = indexw[wbegin] - 1;
//perform limited bisection (if necessary) to get approximate
//eigenvalues to the precision needed.
		    Rlarrb(in, &d[ibegin], &work[indlld + ibegin - 1], p,
			   q, rtol1, rtol2, offset, &work[wbegin], &wgap[wbegin], &werr[wbegin], &work[indwrk], &iwork[iindwk], pivmin, spdiam, in, &iinfo);
		    if (iinfo != 0) {
			*info = -1;
			return;
		    }
//We also recompute the extremal gaps. W holds all eigenvalues
//of the unshifted matrix and must be used for computation
//of WGAP, the entries of WORK might stem from RRRs with
//different shifts. The gaps from WBEGIN-1+OLDFST to
//WBEGIN-1+OLDLST are correctly computed in DLARRB.
//However, we only allow the gaps to become greater since
//this is what should happen when we decrease WERR
		    if (oldfst > 1) {
			mtemp1 = wgap[wbegin + oldfst - 2];
			mtemp2 = w[wbegin + oldfst - 1] - werr[wbegin + oldfst - 1] - w[wbegin + oldfst - 2] - werr[wbegin + oldfst - 2];
			wgap[wbegin + oldfst - 2] = max(mtemp1, mtemp2);
		    }
		    if (wbegin + oldlst - 1 < wend) {
			mtemp1 = wgap[wbegin + oldlst - 1];
			mtemp2 = w[wbegin + oldlst] - werr[wbegin + oldlst] - w[wbegin + oldlst - 1] - werr[wbegin + oldlst - 1];
			wgap[wbegin + oldlst - 1] = max(mtemp1, mtemp2);
		    }
//Each time the eigenvalues in WORK get refined, we store
//the newly found approximation with all shifts applied in W
		    for (j = oldfst; j <= oldlst; j++) {
			w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
		    }
		}
//Process the current node.
		newfst = oldfst;
		for (j = oldfst; j <= oldlst; j++) {
		    if (j == oldlst) {
//we are at the right end of the cluster, this is also the
//boundary of the child cluster
			newlst = j;
		    } else if (wgap[wbegin + j - 1] >= minrgp * abs(work[wbegin + j - 1])) {
//the right relative gap is big enough, the child cluster
//(NEWFST,..,NEWLST) is well separated from the following
			newlst = j;
		    } else {
//inside a child cluster, the relative gap is not
//big enough.
			goto L140;
		    }
//Compute size of child cluster found
		    newsiz = newlst - newfst + 1;
//NEWFTT is the place in Z where the new RRR or the computed
//eigenvector is to be stored
		    if (dol == 1 && dou == m) {
//Store representation at location of the leftmost evalue
//of the cluster
			newftt = wbegin + newfst - 1;
		    } else {
			if (wbegin + newfst - 1 < dol) {
//Store representation at the left end of Z array
			    newftt = dol - 1;
			} else if (wbegin + newfst - 1 > dou) {
//Store representation at the right end of Z array
			    newftt = dou;
			} else {
			    newftt = wbegin + newfst - 1;
			}
		    }
		    if (newsiz > 1) {
//Current child is not a singleton but a cluster.
//Compute and store new representation of child.
//Compute left and right cluster gap.
//LGAP and RGAP are not computed from WORK because
//the eigenvalue approximations may stem from RRRs
//different shifts. However, W hold all eigenvalues
//of the unshifted matrix. Still, the entries in WGAP
//have to be computed from WORK since the entries
//in W might be of the same order so that gaps are not
//exhibited correctly for very close eigenvalues.
			if (newfst == 1) {
			    mtemp1 = Zero, mtemp2 = w[wbegin] - werr[wbegin] - vl;
			    lgap = max(mtemp1, mtemp2);
			} else {
			    lgap = wgap[wbegin + newfst - 2];
			}
			rgap = wgap[wbegin + newlst - 1];
//Compute left- and rightmost eigenvalue of child
//to high precision in order to shift as close
//as possible and obtain as large relative gaps
//as possible
			for (k = 0; k < 2; k++) {
			    if (k == 1) {
				p = indexw[wbegin - 1 + newfst];
			    } else {
				p = indexw[wbegin - 1 + newlst];
			    }
			    offset = indexw[wbegin] - 1;
			    Rlarrb(in, &d[ibegin], &work[indlld + ibegin - 1], p, p, rqtol, rqtol, offset, &work[wbegin], &wgap[wbegin], &werr[wbegin]
				   , &work[indwrk], &iwork[iindwk], pivmin, spdiam, in, &iinfo);
			}
			if (wbegin + newlst - 1 < dol || wbegin + newfst - 1 > dou) {
//if the cluster contains no desired eigenvalues
//skip the computation of that branch of the rep. tree
//We could skip before the refinement of the extremal
//eigenvalues of the child, but then the representation
//tree could be different from the one when nothing is
//skipped. For this reason we skip at this place.
			    idone = idone + newlst - newfst + 1;
			    goto L139;
			}
//Compute RRR of child cluster.
//Note that the new RRR is stored in Z
//DLARRF needs LWORK = 2*N
			Rlarrf(in, &d[ibegin], &l[ibegin], &work[indld +
								 ibegin - 1], newfst, newlst, &work[wbegin],
			       &wgap[wbegin], &werr[wbegin], spdiam, lgap, rgap, pivmin, &tau, &work[indin1], &work[indin2], &work[indwrk], &iinfo);
//In the complex case, DLARRF cannot write
//the new RRR directly into Z and needs an intermediate
//workspace
			for (k = 0; k < in - 1; k++) {
			    z[ibegin + k - 1 + newftt * ldz] = work[indin1 + k - 1];
			    z[ibegin + k - 1 + (newftt + 1) * ldz] = work[indin2 + k - 1];
			}
			z[iend + newftt * ldz] = work[indin1 + in - 1];
			if (iinfo == 0) {
//a new RRR for the cluster was found by DLARRF
//update shift and store it
			    ssigma = sigma + tau;
			    z[iend + (newftt + 1) * ldz] = ssigma;
//WORK() are the midpoints and WERR() the semi-width
//Note that the entries in W are unchanged.
			    for (k = newfst; k <= newlst; k++) {
				fudge = eps * Three * abs(work[wbegin + k - 1]);
				work[wbegin + k - 1] = work[wbegin + k - 1] - tau;
				fudge = fudge + eps * Four * abs(work[wbegin + k - 1]);
//Fudge errors
				werr[wbegin + k - 1] = werr[wbegin + k - 1] + fudge;
//Gaps are not fudged. Provided that WERR is small
//when eigenvalues are close, a zero gap indicates
//that a new representation is needed for resolving
//the cluster. A fudge could lead to a wrong decision
//of judging eigenvalues 'separated' which in
//reality are not. This could have a negative impact
//on the orthogonality of the computed eigenvectors.
			    }
			    ++nclus;
			    k = newcls + (nclus * 2);
			    iwork[k - 1] = newfst;
			    iwork[k] = newlst;
			} else {
			    *info = -2;
			    return;
			}
		    } else {
//Compute eigenvector of singleton
			iter = 0;
			tol = log((double) in) * Four * eps;
			k = newfst;
			windex = wbegin + k - 1;
			windmn = max(windex - 1, (INTEGER) 1);
			windpl = min(windex + 1, m);
			lambda = work[windex];
			++done;
//Check if eigenvector computation is to be skipped
			if (windex < dol || windex > dou) {
			    eskip = MTRUE;
			    goto L125;
			} else {
			    eskip = MFALSE;
			}
			left = work[windex] - werr[windex];
			right = work[windex] + werr[windex];
			indeig = indexw[windex];
//Note that since we compute the eigenpairs for a child,
//all eigenvalue approximations are w.r.t the same shift.
//In this case, the entries in WORK should be used for
//computing the gaps since they exhibit even very small
//differences in the eigenvalues, as opposed to the
//entries in W which might "look" the same.
			if (k == 1) {
//In the case RANGE='I' and with not much initial
//accuracy in LAMBDA and VL, the formula
//LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA )
//can lead to an overestimation of the left gap and
//thus to inadequately early RQI 'convergence'.
//Prevent this by forcing a small left gap.
			    mtemp1 = abs(left), mtemp2 = abs(right);
			    lgap = eps * max(mtemp1, mtemp2);
			} else {
			    lgap = wgap[windmn];
			}
			if (k == im) {
//In the case RANGE='I' and with not much initial
//accuracy in LAMBDA and VU, the formula
//can lead to an overestimation of the right gap and
//thus to inadequately early RQI 'convergence'.
//Prevent this by forcing a small right gap.
			    mtemp1 = abs(left), mtemp2 = abs(right);
			    rgap = eps * max(mtemp1, mtemp2);
			} else {
			    rgap = wgap[windex];
			}
			gap = min(lgap, rgap);
			if (k == 1 || k == im) {
//The eigenvector support can become wrong
//because significant entries could be cut off due to a
//large GAPTOL parameter in LAR1V. Prevent this.
			    gaptol = Zero;
			} else {
			    gaptol = gap * eps;
			}
			isupmn = in;
			isupmx = 1;
//Update WGAP so that it holds the minimum gap
//to the left or the right. This is crucial in the
//case where bisection is used to ensure that the
//eigenvalue is refined up to the required precision.
//The correct value is restored afterwards.
			savgap = wgap[windex];
			wgap[windex] = gap;
//We want to use the Rayleigh Quotient Correction
//as often as possible since it converges quadratically
//when we are close enough to the desired eigenvalue.
//However, the Rayleigh Quotient can have the wrong sign
//and lead us away from the desired eigenvalue. In this
//case, the best we can do is to use bisection.
			usedbs = MFALSE;
			usedrq = MFALSE;
//Bisection is initially turned off unless it is forced
			needbs = !tryrqc;
		      L120:
//Check if bisection should be used to refine eigenvalue
			if (needbs) {
//Take the bisection as new iterate
			    usedbs = MTRUE;
			    itmp1 = iwork[iindr + windex];
			    offset = indexw[wbegin] - 1;
			    Rlarrb(in, &d[ibegin], &work[indlld + ibegin
							 - 1], indeig, indeig, Zero, eps * Two, offset, &work[wbegin],
				   &wgap[wbegin], &werr[wbegin], &work[indwrk], &iwork[iindwk], pivmin, spdiam, itmp1, &iinfo);
			    if (iinfo != 0) {
				*info = -3;
				return;
			    }
			    lambda = work[windex];
//Reset twist index from inaccurate LAMBDA to
//force computation of true MINGMA
			    iwork[iindr + windex] = 0;
			}
//Given LAMBDA, compute the eigenvector.
			Clar1v(in, 1, in, lambda, &d[ibegin], &l[ibegin], &work[indld + ibegin - 1], &work[indlld + ibegin - 1],
			       pivmin, gaptol, &z[ibegin + windex * ldz], !usedbs, &negcnt, &ztz, &mingma, &iwork[iindr + windex],
			       &isuppz[(windex * 2) - 1], &nrminv, &resid, &rqcorr, &work[indwrk]);
			if (iter == 0) {
			    bstres = resid;
			    bstw = lambda;
			} else if (resid < bstres) {
			    bstres = resid;
			    bstw = lambda;
			}
			isupmn = min(isupmn, isuppz[(windex * 2) - 1]);
			isupmx = max(isupmx, isuppz[windex * 2]);
			iter++;
//sin alpha <= |resid|/gap
//Note that both the residual and the gap are
//proportional to the matrix, so ||T|| doesn't play
//a role in the quotient
//Convergence test for Rayleigh-Quotient iteration
//(omitted when Bisection has been used)

			if (resid > tol * gap && abs(rqcorr) > rqtol * abs(lambda) && !usedbs) {
//We need to check that the RQCORR update doesn't
//move the eigenvalue away from the desired one and
//towards a neighbor. -> protection with bisection
			    if (indeig <= negcnt) {
//The wanted eigenvalue lies to the left
				sgndef = -One;
			    } else {
//The wanted eigenvalue lies to the right
				sgndef = One;
			    }
//We only use the RQCORR if it improves the
//the iterate reasonably.
			    if (rqcorr * sgndef >= Zero && lambda + rqcorr <= right && lambda + rqcorr >= left) {
				usedrq = MTRUE;
//Store new midpoint of bisection interval in WORK
				if (sgndef == One) {
//The current LAMBDA is on the left of the true
//eigenvalue
				    left = lambda;
//We prefer to assume that the error estimate
//is correct. We could make the interval not
//as a bracket but to be modified if the RQCORR
//chooses to. In this case, the RIGHT side should
//be modified as follows:
// RIGHT = MAX(RIGHT, LAMBDA + RQCORR)
				} else {
//The current LAMBDA is on the right of the true
//eigenvalue
				    right = lambda;
//ee comment about assuming the error estimate is
//orrect above.
//LEFT = MIN(LEFT, LAMBDA + RQCORR)
				}
				work[windex] = (right + left) * Half;
//Take RQCORR since it has the correct sign and
//improves the iterate reasonably
				lambda = lambda + rqcorr;
//Update width of error interval
				werr[windex] = (right - left) * Half;
			    } else {
				needbs = MTRUE;
			    }
			    if (right - left < rqtol * abs(lambda)) {
//The eigenvalue is computed to bisection accuracy
//compute eigenvector and stop
				usedbs = MTRUE;
				goto L120;
			    } else if (iter < 10) {
				goto L120;
			    } else if (iter == 10) {
				needbs = MTRUE;
				goto L120;
			    } else {
				*info = 5;
				return;
			    }
			} else {
			    stp2ii = MFALSE;
			    if (usedrq && usedbs && bstres <= resid) {
				lambda = bstw;
				stp2ii = MTRUE;
			    }
			    if (stp2ii) {
//improve error angle by second step
				Clar1v(in, 1, in, lambda, &d[ibegin], &l[ibegin], &work[indld + ibegin - 1],
				       &work[indlld + ibegin - 1], pivmin, gaptol, &z[ibegin + windex
										      * ldz], !usedbs, &negcnt, &ztz, &mingma,
				       &iwork[iindr + windex], &isuppz[(windex * 2) - 1], &nrminv, &resid, &rqcorr, &work[indwrk]);
			    }
			    work[windex] = lambda;
			}
//Compute FP-vector support w.r.t. whole matrix
			isuppz[(windex * 2) - 1] = isuppz[(windex * 2) - 1] + oldien;
			isuppz[windex * 2] = isuppz[windex * 2] + oldien;
			zfrom = isuppz[(windex * 2) - 1];
			zto = isuppz[windex * 2];
			isupmn = isupmn + oldien;
			isupmx = isupmx + oldien;
//Ensure vector is ok if support in the RQI has changed
			if (isupmn < zfrom) {
			    for (ii = isupmn; ii <= zfrom - 1; ii++) {
				z[ii + windex * ldz] = Zero;
			    }
			}
			if (isupmx > zto) {
			    for (ii = zto + 1; ii <= isupmx; ii++) {
				z[ii + windex * ldz] = Zero;
			    }
			}
			CRscal(zto - zfrom + 1, nrminv, &z[zfrom + windex * ldz], 1);
		      L125:
//Update W
			w[windex] = lambda + sigma;
//Recompute the gaps on the left and right
//But only allow them to become larger and not
//smaller (which can only happen through "bad"
//cancellation and doesn't reflect the theory
//where the initial gaps are underestimated due
//to WERR being too crude.)
			if (!eskip) {
			    if (k > 1) {
				mtemp1 = wgap[windmn], mtemp2 = w[windex] - werr[windex] - w[windmn] - werr[windmn];
				wgap[windmn] = max(mtemp1, mtemp2);
			    }
			    if (windex < wend) {
				mtemp1 = savgap, mtemp2 = w[windpl] - werr[windpl] - w[windex] - werr[windex];
				wgap[windex] = max(mtemp1, mtemp2);
			    }
			}
			idone++;
		    }
//here ends the code for the current child
		  L139:
//Proceed to any remaining child nodes
		    newfst = j + 1;
		  L140:
		    ;
		}
	    }
	    ++ndepth;
	    goto L40;
	}
	ibegin = iend + 1;
	wbegin = wend + 1;
      L170:
	;
    }
    return;
}
