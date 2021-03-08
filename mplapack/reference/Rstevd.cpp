/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rstevd.cpp,v 1.10 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rstevd(const char *jobz, INTEGER n, REAL * d, REAL * e, REAL * z, INTEGER ldz, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    REAL eps, rmin, rmax, tnrm;
    REAL sigma = 0.0;
    INTEGER lwmin;
    INTEGER wantz;
    INTEGER iscale;
    REAL safmin;
    REAL bignum;
    INTEGER liwmin;
    REAL smlnum;
    INTEGER lquery;
    REAL One = 1.0, Zero = 0.0;

    wantz = Mlsame(jobz, "V");
    lquery = 0;
    if (lwork == -1 || liwork == -1)
	lquery = 1;

    *info = 0;
    liwmin = 1;
    lwmin = 1;
    if (n > 1 && wantz) {
	lwmin = (n * 2) + 1 + n * n;
	liwmin = n * 5 + 3;
    }

    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -6;
    }
    if (*info == 0) {
	work[1] = lwmin;
	iwork[1] = liwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -8;
	} else if (liwork < liwmin && !lquery) {
	    *info = -10;
	}
    }
    if (*info != 0) {
	Mxerbla("Rstevd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
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
    iscale = 0;
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
    }
//For eigenvalues only, call DSTERF.  For eigenvalues and
//eigenvectors, call DSTEDC.
    if (!wantz) {
	Rsterf(n, &d[0], &e[0], info);
    } else {
	Rstedc("I", n, &d[0], &e[0], &z[0], ldz, &work[0], lwork, &iwork[1], liwork, info);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
    if (iscale == 1) {
	Rscal(n, One / sigma, &d[0], 1);
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
