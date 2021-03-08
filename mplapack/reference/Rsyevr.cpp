/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsyevr.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsyevr(const char *jobz, const char *range, const char *uplo, INTEGER n,
	    REAL * A, INTEGER lda, REAL vl, REAL vu, INTEGER il, INTEGER iu, REAL abstol, INTEGER * m, REAL * w,
	    REAL * z, INTEGER ldz, INTEGER * isuppz, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, j, nb, jj;
    REAL eps, vll, vuu, tmp1;
    INTEGER indd, inde;
    REAL anrm;
    INTEGER imax;
    REAL rmin, rmax;
    INTEGER inddd, indee;
    REAL sigma = 0.0;
    INTEGER iinfo;
    char order;
    INTEGER indwk;
    INTEGER lwmin;
    LOGICAL lower, wantz;
    LOGICAL alleig, indeig;
    INTEGER iscale, ieeeok, indibl, indifl;
    LOGICAL valeig;
    REAL safmin;
    REAL abstll, bignum;
    INTEGER indtau, indisp;
    INTEGER indiwo, indwkn;
    INTEGER liwmin;
    LOGICAL tryrac;
    INTEGER llwrkn, llwork, nsplit;
    REAL smlnum;
    INTEGER lwkopt;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    ieeeok = iMlaenv(10, "DSYEVR", "N", 1, 2, 3, 4);
    lower = Mlsame(uplo, "L");
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    lquery = lwork == -1 || liwork == -1;
    lwmin = max((INTEGER) 1, n * 26);
    liwmin = max((INTEGER) 1, n * 10);
    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(alleig || valeig || indeig)) {
	*info = -2;
    } else if (!(lower || Mlsame(uplo, "U"))) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else {
	if (valeig) {
	    if (n > 0 && vu <= vl) {
		*info = -8;
	    }
	} else if (indeig) {
	    if (il < 1 || il > max((INTEGER) 1, n)) {
		*info = -9;
	    } else if (iu < min(n, il) || iu > n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
	if (ldz < 1 || (wantz && ldz < n)) {
	    *info = -15;
	} else if (lwork < lwmin && !lquery) {
	    *info = -18;
	} else if (liwork < liwmin && !lquery) {
	    *info = -20;
	}
    }
    if (*info == 0) {
	nb = iMlaenv(1, "Rsytrd", uplo, n, -1, -1, -1);
	nb = max(nb, iMlaenv(1, "Rormtr", uplo, n, -1, -1, -1));
	lwkopt = max((nb + 1) * n, lwmin);
	work[1] = lwkopt;
	iwork[1] = liwmin;
    }
    if (*info != 0) {
	Mxerbla("Rsyevr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    m = 0;
    if (n == 0) {
	work[1] = One;
	return;
    }
    if (n == 1) {
	work[1] = 7.;
	if (alleig || indeig) {
	    (*m) = 1;
	    w[1] = A[lda + 1];
	} else {
	    if (vl < A[lda + 1] && vu >= A[lda + 1]) {
		(*m) = 1;
		w[1] = A[lda + 1];
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
    abstll = abstol;
    vll = vl;
    vuu = vu;
    anrm = Rlansy("M", uplo, n, &A[0], lda, &work[0]);
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    for (j = 0; j < n; j++) {
		Rscal(n - j + 1, sigma, &A[j + j * lda], 1);
	    }
	} else {
	    for (j = 0; j < n; j++) {
		Rscal(j, sigma, &A[j * lda + 1], 1);
	    }
	}
	if (abstol > Zero) {
	    abstll = abstol * sigma;
	}
	if (valeig) {
	    vll = vl * sigma;
	    vuu = vu * sigma;
	}
    }
//Initialize indices into workspaces.  Note: The IWORK indices are
//used only if DSTERF or DSTEMR fail.
//WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
//elementary reflectors used in DSYTRD.
    indtau = 1;
//WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
    indd = indtau + n;
//WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
//tridiagonal matrix from DSYTRD.
    inde = indd + n;
//WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
//-written by DSTEMR (the DSTERF path copies the diagonal to W).
    inddd = inde + n;
//WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
//-written while computing the eigenvalues in DSTERF and DSTEMR.
    indee = inddd + n;
//INDWK is the starting offset of the left-over workspace, and
//LLWORK is the remaining workspace size.
    indwk = indee + n;
    llwork = lwork - indwk + 1;
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
//Call DSYTRD to reduce symmetric matrix to tridiagonal form.
    Rsytrd(uplo, n, &A[0], lda, &work[indd], &work[inde], &work[indtau], &work[indwk], llwork, &iinfo);
//If all eigenvalues are desired
//then call DSTERF or DSTEMR and DORMTR.
    if ((alleig || (indeig && il == 1 && iu == n)) && ieeeok == 1) {
	if (!wantz) {
	    Rcopy(n, &work[indd], 1, &w[1], 1);
	    Rcopy(n - 1, &work[inde], 1, &work[indee], 1);
	    Rsterf(n, &w[1], &work[indee], info);
	} else {
	    Rcopy(n - 1, &work[inde], 1, &work[indee], 1);
	    Rcopy(n, &work[indd], 1, &work[inddd], 1);
	    if (abstol <= n * Zero * eps) {
		tryrac = MTRUE;
	    } else {
		tryrac = MFALSE;
	    }
	    Rstemr(jobz, "A", n, &work[inddd], &work[indee], vl, vu, il, iu, m, &w[1], &z[0], ldz, n, &isuppz[1], &tryrac, &work[indwk], lwork, &iwork[1], liwork, info);
//Apply orthogonal matrix used in reduction to tridiagonal
//form to eigenvectors returned by DSTEIN.
	    if (wantz && *info == 0) {
		indwkn = inde;
		llwrkn = lwork - indwkn + 1;
		Rormtr("L", uplo, "N", n, (*m), &A[0], lda, &work[indtau], &z[0], ldz, &work[indwkn], llwrkn, &iinfo);
	    }
	}
	if (*info == 0) {
//Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are
//undefined.
	    (*m) = n;
	    goto L30;
	}
	*info = 0;
    }
//Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
//Also call DSTEBZ and DSTEIN if DSTEMR fails.
    if (wantz) {
	order = 'B';
    } else {
	order = 'E';
    }
    Rstebz(range, (const char *) order, n, vll, vuu, il, iu, abstll, &work[indd], &work[inde], m, &nsplit, &w[1], &iwork[indibl],
	   &iwork[indisp], &work[indwk], &iwork[indiwo], info);
    if (wantz) {
	Rstein(n, &work[indd], &work[inde], (*m), &w[1], &iwork[indibl], &iwork[indisp], &z[0], ldz, &work[indwk], &iwork[indiwo], &iwork[indifl], info);
//Apply orthogonal matrix used in reduction to tridiagonal
//form to eigenvectors returned by DSTEIN.
	indwkn = inde;
	llwrkn = lwork - indwkn + 1;
	Rormtr("L", uplo, "N", n, (*m), &A[0], lda, &work[indtau], &z[0], ldz, &work[indwkn], llwrkn, &iinfo);
    }
//  If matrix was scaled, then rescale eigenvalues appropriately.
//Jump here if DSTEMR/DSTEIN succeeded.
  L30:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = (*m);
	} else {
	    imax = *info - 1;
	}
	Rscal(imax, One / sigma, &w[1], 1);
    }
//If eigenvalues are not in order, then sort them, along with
//eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
//It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do
//not return this detailed information to the user.
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
		w[i] = w[j];
		w[j] = tmp1;
		Rswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
	    }
	}
    }
//Set WORK(1) to optimal workspace size.
    work[1] = lwkopt;
    iwork[1] = liwmin;
    return;
}
