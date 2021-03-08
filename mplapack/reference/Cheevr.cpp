/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cheevr.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cheevr(const char *jobz, const char *range, const char *uplo, INTEGER n,
	    COMPLEX * A, INTEGER lda, REAL vl, REAL vu,
	    INTEGER il, INTEGER iu, REAL abstol, INTEGER * m, REAL *
	    w, COMPLEX * z, INTEGER ldz, INTEGER * isuppz, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER lrwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, j, nb, jj;
    REAL eps, vll = 0.0, vuu = 0.0, tmp1, anrm;
    INTEGER imax;
    REAL rmin, rmax;
    LOGICAL test;
    INTEGER itmp1;
    INTEGER indrd, indre;
    REAL sigma = 0.0;
    INTEGER iinfo;
    char order;
    INTEGER indwk;
    INTEGER lwmin;
    LOGICAL lower, wantz;
    LOGICAL alleig, indeig;
    INTEGER iscale, ieeeok, indibl, indrdd, indifl, indree;
    LOGICAL valeig;
    REAL safmin;
    REAL abstll, bignum;
    INTEGER indtau, indisp;
    INTEGER indiwo, indwkn;
    INTEGER indrwk, liwmin;
    LOGICAL tryrac;
    INTEGER lrwmin, llwrkn, llwork, nsplit;
    REAL smlnum;
    LOGICAL lquery;
    INTEGER lwkopt;
    INTEGER llrwork;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    ieeeok = iMlaenv(10, "ZHEEVR", "N", 1, 2, 3, 4);
    lower = Mlsame(uplo, "L");
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;
    lrwmin = max((INTEGER) 1, n * 24);
    liwmin = max((INTEGER) 1, n * 10);
    lwmin = max((INTEGER) 1, n << 1);
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
	}
    }
    if (*info == 0) {
	nb = iMlaenv(1, "ZHETRD", uplo, n, -1, -1, -1);
	nb = max(nb, iMlaenv(1, "ZUNMTR", uplo, n, -1, -1, -1));
	lwkopt = max((nb + 1) * n, lwmin);
	work[1] = lwkopt;
	rwork[1] = lrwmin;
	iwork[1] = liwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -18;
	} else if (lrwork < lrwmin && !lquery) {
	    *info = -20;
	} else if (liwork < liwmin && !lquery) {
	    *info = -22;
	}
    }
    if (*info != 0) {
	Mxerbla("Cheevr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    m = 0;
    if (n == 0) {
	work[1] = 1;
	return;
    }
    if (n == 1) {
	work[1] = 2;
	if (alleig || indeig) {
	    (*m) = 1;
	    w[1] = A[lda + 1].real();
	} else {
	    if (vl < A[lda + 1].real() && vu >= A[lda + 1].real()) {
		(*m) = 1;
		w[1] = A[lda + 1].real();
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
    if (valeig) {
	vll = vl;
	vuu = vu;
    }
    anrm = Clansy("M", uplo, n, &A[0], lda, &rwork[1]);
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
		CRscal(n - j + 1, sigma, &A[j + j * lda], 1);
	    }
	} else {
	    for (j = 0; j < n; j++) {
		CRscal(j, sigma, &A[j * lda + 1], 1);
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
//used only if DSTERF or ZSTEMR fail.
//WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
//elementary reflectors used in ZHETRD.
    indtau = 1;
//INDWK is the starting offset of the remaining complex workspace,
//and LLWORK is the remaining complex workspace size.
    indwk = indtau + n;
    llwork = lwork - indwk + 1;
//RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
//entries.
    indrd = 1;
//RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
//tridiagonal matrix from ZHETRD.
    indre = indrd + n;
//RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
//-written by ZSTEMR (the DSTERF path copies the diagonal to W).
    indrdd = indre + n;
//RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
//-written while computing the eigenvalues in DSTERF and ZSTEMR.
    indree = indrdd + n;
//INDRWK is the starting offset of the left-over real workspace, and
//LLRWORK is the remaining workspace size.
    indrwk = indree + n;
    llrwork = lrwork - indrwk + 1;
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
//Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
    Chetrd(uplo, n, &A[0], lda, &rwork[indrd], &rwork[indre], &work[indtau], &work[indwk], llwork, &iinfo);
//If all eigenvalues are desired
//then call DSTERF or ZSTEMR and ZUNMTR.
    test = MFALSE;
    if (indeig) {
	if (il == 1 && iu == n) {
	    test = MTRUE;
	}
    }
    if ((alleig || test) && ieeeok == 1) {
	if (!wantz) {
	    Rcopy(n, &rwork[indrd], 1, &w[1], 1);
	    Rcopy(n - 1, &rwork[indre], 1, &rwork[indree], 1);
	    Rsterf(n, &w[1], &rwork[indree], info);
	} else {
	    Rcopy(n - 1, &rwork[indre], 1, &rwork[indree], 1);
	    Rcopy(n, &rwork[indrd], 1, &rwork[indrdd], 1);
	    if (abstol <= n * Two * eps) {
		tryrac = MTRUE;
	    } else {
		tryrac = MFALSE;
	    }
	    Cstemr(jobz, "A", n, &rwork[indrdd], &rwork[indree], vl, vu, il, iu, m, &w[1], &z[0], ldz, n, &isuppz[1], &tryrac, &rwork[indrwk], llrwork, &iwork[1], liwork, info);
//Apply unitary matrix used in reduction to tridiagonal
//form to eigenvectors returned by ZSTEIN.
	    if (wantz && *info == 0) {
		indwkn = indwk;
		llwrkn = lwork - indwkn + 1;
		Cunmtr("L", uplo, "N", n, (*m), &A[0], lda, &work[indtau], &z[0], ldz, &work[indwkn], llwrkn, &iinfo);
	    }
	}
	if (*info == 0) {
	    (*m) = n;
	    goto L30;
	}
	*info = 0;
    }
//Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
//Also call DSTEBZ and ZSTEIN if ZSTEMR fails.
    if (wantz) {
	order = 'B';
    } else {
	order = 'E';
    }
    Rstebz(range, (const char *) order, n, vll, vuu, il, iu, abstll, &rwork[indrd], &rwork[indre], m, &nsplit, &w[1],
	   &iwork[indibl], &iwork[indisp], &rwork[indrwk], &iwork[indiwo], info);
    if (wantz) {
	Cstein(n, &rwork[indrd], &rwork[indre], (*m), &w[1], &iwork[indibl], &iwork[indisp], &z[0], ldz, &rwork[indrwk], &iwork[indiwo], &iwork[indifl], info);
//Apply unitary matrix used in reduction to tridiagonal
//form to eigenvectors returned by ZSTEIN.
	indwkn = indwk;
	llwrkn = lwork - indwkn + 1;
	Cunmtr("L", uplo, "N", n, (*m), &A[0], lda, &work[indtau], &z[0], ldz, &work[indwkn], llwrkn, &iinfo);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
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
		itmp1 = iwork[indibl + i - 1];
		w[i] = w[j];
		iwork[indibl + i - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		Cswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
	    }
	}
    }
//Set WORK(1) to optimal workspace size.
    work[1] = lwkopt;
    rwork[1] = lrwmin;
    iwork[1] = liwmin;
    return;
}
