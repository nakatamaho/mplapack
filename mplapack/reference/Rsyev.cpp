/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsyev.cpp,v 1.11 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsyev(const char *jobz, const char *uplo, INTEGER n, REAL * A, INTEGER lda, REAL * w, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER wantz, lower, lquery, nb, lwkopt, iscale, imax;
    INTEGER inde, indtau, indwrk, llwork, iinfo;
    REAL safmin, eps, smlnum, bignum, rmin, rmax;
    REAL sigma = 0.0, anrm;
    REAL rtmp;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    lower = Mlsame(uplo, "L");
    lquery = lwork == -1;
    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(lower || Mlsame(uplo, "U"))) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    }
    if (*info == 0) {
	nb = iMlaenv(1, "Rsytrd", uplo, n, -1, -1, -1);
	lwkopt = max((INTEGER) 1, (nb + 2) * n);
	work[0] = lwkopt;
	if (lwork < max((INTEGER) 1, 3 * n - 1) && !lquery) {
	    *info = -8;
	}
    }
    if (*info != 0) {
	Mxerbla("Rsyev ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
	w[0] = A[0];
	work[0] = Two;
	if (wantz) {
	    A[0] = One;
	}
	return;
    }
//Get machine constants.
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
    smlnum = safmin / eps;
    bignum = One / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

//Scale matrix to allowable range, if necessary.
    anrm = Rlansy("M", uplo, n, A, lda, work);
    iscale = 0;
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	Rlascl(uplo, 0, 0, One, sigma, n, n, A, lda, info);
    }
//Call Rsytrd to reduce symmetric matrix to tridiagonal form.
    inde = 1;
    indtau = inde + n;
    indwrk = indtau + n;
    llwork = lwork - indwrk + 1;
    Rsytrd(uplo, n, A, lda, w, &work[inde - 1], &work[indtau - 1], &work[indwrk - 1], llwork, &iinfo);
//For eigenvalues only, call Rsterf.  For eigenvectors, first call
//Rorgtr to generate the orthogonal matrix, then call Rsteqr.
    if (!wantz) {
	Rsterf(n, w, &work[inde - 1], info);
    } else {
	Rorgtr(uplo, n, A, lda, &work[indtau - 1], &work[indwrk - 1], llwork, &iinfo);
	Rsteqr(jobz, n, w, &work[inde - 1], A, lda, &work[indtau - 1], info);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
    if (iscale == 1) {
	if (*info == 0) {
	    imax = n;
	} else {
	    imax = *info - 1;
	}
	rtmp = One / sigma;
	Rscal(imax, rtmp, w, 1);
    }
//Set WORK(0) to optimal workspace size.
    work[0] = lwkopt;
    return;
}
