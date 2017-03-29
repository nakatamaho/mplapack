/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rstevr.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rstevr(const char *jobz, const char *range, INTEGER n, REAL *
	    d, REAL * e, REAL vl, REAL vu, INTEGER il, INTEGER iu, REAL abstol, INTEGER * m, REAL * w, REAL * z, INTEGER ldz,
	    INTEGER * isuppz, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, j, jj;
    REAL eps, vll, vuu, tmp1;
    INTEGER imax;
    REAL rmin, rmax;
    LOGICAL test;
    REAL tnrm;
    INTEGER itmp1;
    REAL sigma = 0.0;
    char order;
    INTEGER lwmin, indiwo;
    LOGICAL wantz;
    LOGICAL alleig, indeig;
    INTEGER iscale, ieeeok, indibl, indifl;
    LOGICAL valeig;
    REAL safmin;
    REAL bignum;
    INTEGER indisp;
    INTEGER liwmin;
    LOGICAL tryrac;
    INTEGER nsplit;
    REAL smlnum;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    ieeeok = iMlaenv(10, "Rstevr", "N", 1, 2, 3, 4);
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");

    lquery = lwork == -1 || liwork == -1;
    lwmin = max((INTEGER) 1, n * 20);
    liwmin = max((INTEGER) 1, n * 10);
    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(alleig || valeig || indeig)) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else {
	if (valeig) {
	    if (n > 0 && vu <= vl) {
		*info = -7;
	    }
	} else if (indeig) {
	    if (il < 1 || il > max((INTEGER) 1, n)) {
		*info = -8;
	    } else if (iu < min(n, il) || iu > n) {
		*info = -9;
	    }
	}
    }
    if (*info == 0) {
	if (ldz < 1 || (wantz && ldz < n)) {
	    *info = -14;
	}
    }
    if (*info == 0) {
	work[1] = lwmin;
	iwork[1] = liwmin;
	if (lwork < (lwmin && !lquery)) {
	    *info = -17;
	} else if (liwork < (liwmin && !lquery)) {
	    *info = -19;
	}
    }
    if (*info != 0) {
	Mxerbla("Rstevr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    (*m) = 0;
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (alleig || indeig) {
	    (*m) = 1;
	    w[1] = d[1];
	} else {
	    if (vl < d[1] && vu >= d[1]) {
		(*m) = 1;
		w[1] = d[1];
	    }
	}
	if (wantz) {
	    z[ldz + 1] = One;
	}
	return;
    }
//Get machine constants.
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
    smlnum = safmin / eps;
    bignum = One / smlnum;
    rmin = sqrt(smlnum);
    mtemp1 = sqrt(bignum), mtemp2 = One / sqrt(sqrt(safmin));
    rmax = min(mtemp1, mtemp2);
//Scale matrix to allowable range, if necessary.
    iscale = 0;
    vll = vl;
    vuu = vu;
    tnrm = Rlanst("M", n, &d[0], &e[0]);
    if (tnrm > Zero && tnrm < rmin) {
	iscale = 1;
	sigma = rmin / tnrm;
    } else if (tnrm > rmax) {
	iscale = 1;
	sigma = rmax / tnrm;
    }
    if (iscale == 1) {
	Rscal(n, sigma, &d[0], 1);
	Rscal(n - 1, sigma, &e[0], 1);
	if (valeig) {
	    vll = vl * sigma;
	    vuu = vu * sigma;
	}
    }
//Initialize indices into workspaces.  Note: These indices are used only
//if DSTERF or DSTEMR fail.
//IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
//stores the block indices of each of the M<=N eigenvalues.
    indibl = 0;
//IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
//stores the starting and finishing indices of each block.
    indisp = indibl + n;
//IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
//that corresponding to eigenvectors that fail to converge in
//DSTEIN.  This information is discarded; if any fail, the driver
//returns INFO > Zero
    indifl = indisp + n;
//INDIWO is the offset of the remaining INTEGER workspace.
    indiwo = indisp + n;
//If all eigenvalues are desired, then
//call DSTERF or DSTEMR.  If this fails for some eigenvalue, then
//try DSTEBZ.
    test = MFALSE;
    if (indeig) {
	if (il == 1 && iu == n) {
	    test = MTRUE;
	}
    }
    if ((alleig || test) && ieeeok == 1) {
	Rcopy(n - 1, &e[0], 1, &work[0], 1);
	if (!wantz) {
	    Rcopy(n, &d[0], 1, &w[0], 1);
	    Rsterf(n, &w[0], &work[0], info);
	} else {
	    Rcopy(n, &d[0], 1, &work[n + 1], 1);
	    if (abstol <= n * Two * eps) {
		tryrac = MTRUE;
	    } else {
		tryrac = MFALSE;
	    }
	    Rstemr(jobz, "A", n, &work[n + 1], &work[0], vl, vu, il, iu, m, &w[0], &z[0], ldz, n, &isuppz[1], &tryrac,
		   &work[(n << 1) + 1], lwork - (n << 1), &iwork[1], liwork, info);
	}
	if (*info == 0) {
	    (*m) = n;
	    goto L10;
	}
	*info = 0;
    }
//Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
    if (wantz) {
	order = 'B';
    } else {
	order = 'E';
    }
    Rstebz(range, (const char *) order, n, vll, vuu, il, iu, abstol, &d[0], &e[0], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[0], &iwork[indiwo], info);
    if (wantz) {
	Rstein(n, &d[0], &e[0], (*m), &w[1], &iwork[indibl], &iwork[indisp], &z[0], ldz, &work[0], &iwork[indiwo], &iwork[indifl], info);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
  L10:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = (*m);
	} else {
	    imax = *info - 1;
	}
	Rscal(imax, One / sigma, &w[1], 1);
    }
//If eigenvalues are not in order, then sort them, along with
//eigenvectors.
    if (wantz) {
	for (j = 0; j < (*m) - 1; j++) {
	    i = 0;
	    tmp1 = w[j];
	    for (jj = j + 1; jj <= (*m); jj++) {
		if (w[jj] < tmp1) {
		    i = jj;
		    tmp1 = w[jj];
		}
	    }
	    if (i != 0) {
		itmp1 = iwork[i];
		w[i] = w[j];
		iwork[i] = iwork[j];
		w[j] = tmp1;
		iwork[j] = itmp1;
		Rswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
	    }
	}
    }
//Causes problems with tests 19 & 20:
//IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1002 + .002
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
