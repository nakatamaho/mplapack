/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chbev.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chbev(const char *jobz, const char *uplo, INTEGER n, INTEGER kd, COMPLEX * AB, INTEGER ldab, REAL * w, COMPLEX * z, INTEGER ldz, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    REAL eps;
    INTEGER inde;
    REAL anrm;
    INTEGER imax;
    REAL rmin, rmax;
    REAL sigma = 0.0;
    INTEGER iinfo;
    INTEGER lower, wantz;
    INTEGER iscale;
    REAL safmin;
    REAL bignum;
    INTEGER indrwk;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    lower = Mlsame(uplo, "L");
    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(lower || Mlsame(uplo, "U"))) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (kd < 0) {
	*info = -4;
    } else if (ldab < kd + 1) {
	*info = -6;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -9;
    }
    if (*info != 0) {
	Mxerbla("Chbev ", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (lower) {
	    w[1] = AB[ldab + 1].real();
	} else {
	    w[1] = AB[kd + 1 + ldab].real();
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
    rmax = sqrt(bignum);
//Scale matrix to allowable range, if necessary.
    anrm = Clanhb("M", uplo, n, kd, AB, ldab, rwork);
    iscale = 0;
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    Clascl("B", kd, kd, One, sigma, n, n, AB, ldab, info);
	} else {
	    Clascl("Q", kd, kd, One, sigma, n, n, AB, ldab, info);
	}
    }
//Call ZHBTRD to reduce Hermitian band matrix to tridiagonal form.
    inde = 1;
    Chbtrd(jobz, uplo, n, kd, AB, ldab, w, &rwork[inde], z, ldz, work, &iinfo);
//For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR.
    if (!wantz) {
	Rsterf(n, w, &rwork[inde], info);
    } else {
	indrwk = inde + n;
	Csteqr(jobz, n, w, &rwork[inde], z, ldz, &rwork[indrwk], info);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
    if (iscale == 1) {
	if (*info == 0) {
	    imax = n;
	} else {
	    imax = *info - 1;
	}
	Rscal(imax, One / sigma, w, 1);
    }
    return;
}
