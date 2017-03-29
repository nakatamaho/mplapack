/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cstemr.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cstemr(const char *jobz, const char *range, INTEGER n, REAL * d, REAL * e, REAL vl, REAL vu, INTEGER il, INTEGER iu,
	    INTEGER * m, REAL * w, COMPLEX * z, INTEGER ldz, INTEGER nzc, INTEGER * isuppz, LOGICAL * tryrac, REAL * work,
	    INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, j;
    REAL r1, r2;
    INTEGER jj;
    REAL cs;
    INTEGER in;
    REAL sn, wl, wu;
    INTEGER iil, iiu;
    REAL eps, tmp;
    INTEGER indd, iend, jblk, wend;
    REAL rmin, rmax;
    INTEGER itmp;
    REAL tnrm;
    INTEGER inde2, itmp2;
    REAL rtol1, rtol2;
    REAL scale;
    INTEGER indgp;
    INTEGER iinfo, iindw, ilast;
    INTEGER lwmin;
    INTEGER wantz;
    INTEGER alleig;
    INTEGER ibegin;
    INTEGER indeig;
    INTEGER iindbl;
    INTEGER valeig;
    INTEGER wbegin;
    REAL safmin;
    REAL bignum;
    INTEGER inderr, iindwk, indgrs, offset;
    REAL thresh;
    INTEGER iinspl, indwrk, ifirst, liwmin, nzcmin;
    REAL pivmin;
    INTEGER nsplit;
    REAL smlnum;
    INTEGER lquery, zquery;
    REAL Zero = 0.0, One = 1.0, MINRGP = 0.001, FTOL = 0.005, Four = 4.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    lquery = lwork == -1 || liwork == -1;
    zquery = nzc == -1;
    *tryrac = *info != 0;
//DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
//In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
//Furthermore, ZLARRV needs WORK of size 12*N, IWORK of size 7*N.
    if (wantz) {
	lwmin = n * 18;
	liwmin = n * 10;
    } else {
//need less workspace if only the eigenvalues are wanted
	lwmin = n * 12;
	liwmin = n * 8;
    }
    wl = Zero;
    wu = Zero;
    iil = 0;
    iiu = 0;
    if (valeig) {
//We do not reference VL, VU in the cases RANGE = 'I','A'
//The INTEGERerval (WL, WU] contains all the wanted eigenvalues.
//It is either given by the user or computed in DLARRE.
	wl = vl;
	wu = vu;
    } else if (indeig) {
//We do not reference IL, IU in the cases RANGE = 'V','A'
	iil = il;
	iiu = iu;
    }
    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(alleig || valeig || indeig)) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (valeig && n > 0 && wu <= wl) {
	*info = -7;
    } else if (indeig && (iil < 1 || iil > n)) {
	*info = -8;
    } else if (indeig && (iiu < iil || iiu > n)) {
	*info = -9;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -13;
    } else if (lwork < lwmin && !lquery) {
	*info = -17;
    } else if (liwork < liwmin && !lquery) {
	*info = -19;
    }
//Get machine constants.
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
    smlnum = safmin / eps;
    bignum = One / smlnum;
    rmin = sqrt(smlnum);
    mtemp1 = sqrt(bignum), mtemp2 = One / sqrt(sqrt(safmin));
    rmax = min(mtemp1, mtemp2);
    if (*info == 0) {
	work[1] = lwmin;
	iwork[1] = liwmin;
	if (wantz && alleig) {
	    nzcmin = n;
	} else if (wantz && valeig) {
	    Rlarrc("T", n, vl, vu, &d[0], &e[0], safmin, &nzcmin, &itmp, &itmp2, info);
	} else if (wantz && indeig) {
	    nzcmin = iiu - iil + 1;
	} else {
//WANTZ .EQ. FALSE.
	    nzcmin = 0;
	}
	if (zquery && *info == 0) {
	    z[ldz + 1] = nzcmin;
	} else if (nzc < nzcmin && !zquery) {
	    *info = -14;
	}
    }
    if (*info != 0) {
	Mxerbla("Cstemr", -(*info));
	return;
    } else if (lquery || zquery) {
	return;
    }
//Handle N = 0, 1, and 2 cases immediately
    m = 0;
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (alleig || indeig) {
	    (*m) = 1;
	    w[1] = d[1];
	} else {
	    if (wl < d[1] && wu >= d[1]) {
		(*m) = 1;
		w[1] = d[1];
	    }
	}
	if (wantz && !zquery) {
	    z[ldz + 1] = One;
	    isuppz[1] = 1;
	    isuppz[2] = 1;
	}
	return;
    }
    if (n == 2) {
	if (!wantz) {
	    Rlae2(d[0], e[0], d[2], &r1, &r2);
	} else if (wantz && !zquery) {
	    Rlaev2(d[0], e[0], d[2], &r1, &r2, &cs, &sn);
	}
	if ((alleig || valeig) && (r2 > wl) && (r2 <= wu || indeig) && iil == 1) {
	    ++(*m);
	    w[*m] = r2;
	    if (wantz && !zquery) {
		z[(*m) * ldz + 1] = -sn;
		z[(*m) * ldz + 2] = cs;
//Note: At most one of SN and CS can be zero.
		if (sn != Zero) {
		    if (cs != Zero) {
			isuppz[((*m) * 2) - 1] = 1;
			isuppz[((*m) * 2) - 1] = 2;
		    } else {
			isuppz[((*m) * 2) - 1] = 1;
			isuppz[((*m) * 2) - 1] = 1;
		    }
		} else {
		    isuppz[((*m) * 2) - 1] = 2;
		    isuppz[((*m) * 2)] = 2;
		}
	    }
	}
	if (alleig || (valeig && r1 > wl && r1 <= wu) || (indeig && iiu == 2)) {
	    ++(*m);
	    w[(*m)] = r1;
	    if (wantz && !zquery) {
		z[(*m) * ldz + 1] = cs;
		z[(*m) * ldz + 2] = sn;
//Note: At most one of SN and CS can be zero.
		if (sn != Zero) {
		    if (cs != Zero) {
			isuppz[((*m) * 2) - 1] = 1;
			isuppz[((*m) * 2) - 1] = 2;
		    } else {
			isuppz[((*m) * 2) - 1] = 1;
			isuppz[((*m) * 2) - 1] = 1;
		    }
		} else {
		    isuppz[((*m) * 2) - 1] = 2;
		    isuppz[(*m) * 2] = 2;
		}
	    }
	}
	return;
    }
//Continue with general N
    indgrs = 1;
    inderr = n * 2 + 1;
    indgp = n * 3 + 1;
    indd = n * 4 + 1;
    inde2 = n * 5 + 1;
    indwrk = n * 6 + 1;
    iinspl = 0;
    iindbl = n + 1;
    iindw = n * 2 + 1;
    iindwk = n * 3 + 1;
//Scale matrix to allowable range, if necessary.
//The allowable range is related to the PIVMIN parameter; see the
//comments in DLARRD.  The preference for scaling small values
//up is heuristic; we expect users' matrices not to be close to the
//RMAX threshold.
    scale = One;
    tnrm = Rlanst("M", n, &d[0], &e[0]);
    if (tnrm > Zero && tnrm < rmin) {
	scale = rmin / tnrm;
    } else if (tnrm > rmax) {
	scale = rmax / tnrm;
    }
    if (scale != One) {
	Rscal(n, scale, &d[0], 1);
	Rscal(n - 1, scale, &e[0], 1);
	tnrm = tnrm * scale;
	if (valeig) {
//If eigenvalues in interval have to be found,
//scale (WL, WU] accordingly
	    wl = wl * scale;
	    wu = wu * scale;
	}
    }
//Compute the desired eigenvalues of the tridiagonal after splitting
//into smaller subblocks if the corresponding off-diagonal elements
//are small
//THRESH is the splitting parameter for DLARRE
//A negative THRESH forces the old splitting criterion based on the
//size of the off-diagonal. A positive THRESH switches to splitting
//which preserves relative accuracy.
    if (*tryrac) {
//Test whether the matrix warrants the more expensive relative approach.
	Rlarrr(n, &d[0], &e[0], &iinfo);
    } else {
//The user does not care about relative accurately eigenvalues
	iinfo = -1;
    }
//Set the splitting criterion
    if (iinfo == 0) {
	thresh = eps;
    } else {
	thresh = -eps;
//relative accuracy is desired but T does not guarantee it
	*tryrac = MFALSE;
    }
    if (*tryrac) {
//Copy original diagonal, needed to guarantee relative accuracy
	Rcopy(n, &d[0], 1, &work[indd], 1);
    }
//Store the squares of the offdiagonal values of T
    for (j = 0; j < n - 1; j++) {
	work[inde2 + j - 1] = e[j] * e[j];
    }
//Set the tolerance parameters for bisection
    if (!wantz) {
//DLARRE computes the eigenvalues to full precision.
	rtol1 = eps * Four;
	rtol2 = eps * Four;
    } else {
//DLARRE computes the eigenvalues to less than full precision.
//ZLARRV will refine the eigenvalue approximations, and we only
//need less accurate initial bisection in DLARRE.
//Note: these settings do only affect the subset case and DLARRE
	rtol1 = sqrt(eps);
	mtemp1 = sqrt(eps) * FTOL, mtemp2 = eps * Four;
	rtol2 = max(mtemp1, mtemp2);
    }
    Rlarre(range, n, &wl, &wu, iil, iiu, &d[0], &e[0], &work[inde2], rtol1, rtol2, thresh, &nsplit, &iwork[iinspl], m, &w[1],
	   &work[inderr], &work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs], &pivmin, &work[indwrk], &iwork[iindwk], &iinfo);
    if (iinfo != 0) {
	*info = abs(iinfo) + 10;
	return;
    }
//Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired
//part of the spectrum. All desired eigenvalues are contained in
//(WL,WU]
    if (wantz) {

//Compute the desired eigenvectors corresponding to the computed
//eigenvalues
	Clarrv(n, wl, wu, &d[0], &e[0], pivmin, &iwork[iinspl], (*m), 1, (*m), MINRGP, rtol1, rtol2, &w[1], &work[inderr],
	       &work[indgp], &iwork[iindbl], &iwork[iindw], &work[indgrs], &z[0], ldz, &isuppz[1], &work[indwrk], &iwork[iindwk], &iinfo);
	if (iinfo != 0) {
	    *info = abs(iinfo) + 20;
	    return;
	}
    } else {
//DLARRE computes eigenvalues of the (shifted) root representation
//ZLARRV returns the eigenvalues of the unshifted matrix.
//However, if the eigenvectors are not desired by the user, we need
//to apply the corresponding shifts from DLARRE to obtain the
//eigenvalues of the original matrix.
	for (j = 0; j < (*m); j++) {
	    itmp = iwork[iindbl + j - 1];
	    w[j] = w[j] + e[iwork[iinspl + itmp - 1]];
	}
    }
    if (*tryrac) {
//Refine computed eigenvalues so that they are relatively accurate
//with respect to the original matrix T.
	ibegin = 1;
	wbegin = 1;
	for (jblk = 0; jblk <= iwork[iindbl + (*m) - 1]; jblk++) {
	    iend = iwork[iinspl + jblk - 1];
	    in = iend - ibegin + 1;
	    wend = wbegin - 1;
//check if any eigenvalues have to be refined in this block
	  L36:
	    if (wend < (*m)) {
		if (iwork[iindbl + wend] == jblk) {
		    ++wend;
		    goto L36;
		}
	    }
	    if (wend < wbegin) {
		ibegin = iend + 1;
		goto L39;
	    }
	    offset = iwork[iindw + wbegin - 1] - 1;
	    ifirst = iwork[iindw + wbegin - 1];
	    ilast = iwork[iindw + wend - 1];
	    rtol2 = eps * Four;
	    Rlarrj(in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 1],
		   ifirst, ilast, rtol2, offset, &w[wbegin], &work[inderr + wbegin - 1], &work[indwrk], &iwork[iindwk], pivmin, tnrm, &iinfo);
	    ibegin = iend + 1;
	    wbegin = wend + 1;
	  L39:
	    ;
	}
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
    if (scale != One) {
	Rscal((*m), One / scale, &w[1], 1);
    }
//If eigenvalues are not in increasing order, then sort them,
//possibly along with eigenvectors.
    if (nsplit > 1) {
	if (!wantz) {
	    Rlasrt("I", (*m), &w[1], &iinfo);
	    if (iinfo != 0) {
		*info = 3;
		return;
	    }
	} else {
	    for (j = 0; j < (*m) - 1; j++) {
		i = 0;
		tmp = w[j];
		for (jj = j + 1; jj <= (*m); jj++) {
		    if (w[jj] < tmp) {
			i = jj;
			tmp = w[jj];
		    }

		}
		if (i != 0) {
		    w[i] = w[j];
		    w[j] = tmp;
		    if (wantz) {
			Cswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
			itmp = isuppz[(i * 2) - 1];
			isuppz[(i * 2) - 1] = isuppz[(j * 2) - 1];
			isuppz[(j * 2) - 1] = itmp;
			itmp = isuppz[i * 2];
			isuppz[i * 2] = isuppz[j * 2];
			isuppz[j * 2] = itmp;
		    }
		}

	    }
	}
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
