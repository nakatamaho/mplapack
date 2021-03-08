/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsbgvx.cpp,v 1.3 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsbgvx(const char *jobz, const char *range, const char *uplo, INTEGER n,
	    INTEGER ka, INTEGER kb, REAL * AB, INTEGER ldab, REAL * BB, INTEGER ldbb, REAL * q, INTEGER ldq, REAL vl, REAL vu,
	    INTEGER il, INTEGER iu, REAL abstol, INTEGER * m, REAL * w, REAL * z, INTEGER ldz, REAL * work, INTEGER * iwork, INTEGER * ifail, INTEGER * info)
{
    INTEGER i, j, jj;
    REAL tmp1;
    INTEGER indd, inde;
    char vect;
    INTEGER test;
    INTEGER itmp1, indee;
    INTEGER iinfo;
    char order;
    INTEGER upper, wantz, alleig, indeig;
    INTEGER indibl;
    INTEGER valeig;
    INTEGER indiwo, indisp;
    INTEGER indwrk;
    INTEGER nsplit;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    upper = Mlsame(uplo, "U");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");

    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(alleig || valeig || indeig)) {
	*info = -2;
    } else if (!(upper || Mlsame(uplo, "L"))) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (ka < 0) {
	*info = -5;
    } else if (kb < 0 || kb > ka) {
	*info = -6;
    } else if (ldab < ka + 1) {
	*info = -8;
    } else if (ldbb < kb + 1) {
	*info = -10;
    } else if (ldq < 1 || (wantz && ldq < n)) {
	*info = -12;
    } else {
	if (valeig) {
	    if (n > 0 && vu <= vl) {
		*info = -14;
	    }
	} else if (indeig) {
	    if (il < 1 || il > max((INTEGER) 1, n)) {
		*info = -15;
	    } else if (iu < min(n, il) || iu > n) {
		*info = -16;
	    }
	}
    }
    if (*info == 0) {
	if (ldz < 1 || (wantz && ldz < n)) {
	    *info = -21;
	}
    }
    if (*info != 0) {
	Mxerbla("Rsbgvx", -(*info));
	return;
    }
//Quick return if possible
    m = 0;
    if (n == 0) {
	return;
    }
//Form a split Cholesky factorization of B.
    Rpbstf(uplo, n, kb, &BB[0], ldbb, info);
    if (*info != 0) {
	*info = n + *info;
	return;
    }
//Transform problem to standard eigenvalue problem.
    Rsbgst(jobz, uplo, n, ka, kb, &AB[0], ldab, &BB[0], ldbb, &q[0], ldq, &work[0], &iinfo);
//Reduce symmetric band matrix to tridiagonal form.
    indd = 1;
    inde = indd + n;
    indwrk = inde + n;
    if (wantz) {
	vect = 'U';
    } else {
	vect = 'N';
    }
    Rsbtrd(&vect, uplo, n, ka, &AB[0], ldab, &work[indd], &work[inde], &q[0], ldq, &work[indwrk], &iinfo);
//If all eigenvalues are desired and ABSTOL is less than or equal
//to zero, then call DSTERF or SSTEQR.  If this fails for some
//eigenvalue, then try DSTEBZ.
    test = MFALSE;
    if (indeig) {
	if (il == 1 && iu == n) {
	    test = MTRUE;
	}
    }
    if ((alleig || test) && abstol <= Zero) {
	Rcopy(n, &work[indd], 1, &w[1], 1);
	indee = indwrk + n * 2;
	Rcopy(n - 1, &work[inde], 1, &work[indee], 1);
	if (!wantz) {
	    Rsterf(n, &w[1], &work[indee], info);
	} else {
	    Rlacpy("A", n, n, &q[0], ldq, &z[0], ldz);
	    Rsteqr(jobz, n, &w[1], &work[indee], &z[0], ldz, &work[indwrk], info);
	    if (*info == 0) {
		for (i = 0; i < n; i++) {
		    ifail[i] = 0;
		}
	    }
	}
	if (*info == 0) {
	    *m = n;
	    goto L30;
	}
	*info = 0;
    }
//Otherwise, call DSTEBZ and, if eigenvectors are desired,
//call DSTEIN.
    if (wantz) {
	order = 'B';
    } else {
	order = 'E';
    }
    indibl = 0;
    indisp = indibl + n;
    indiwo = indisp + n;
    Rstebz(range, &order, n, vl, vu, il, iu, abstol, &work[indd], &work[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk], &iwork[indiwo], info);
    if (wantz) {
	Rstein(n, &work[indd], &work[inde], *m, &w[1], &iwork[indibl], &iwork[indisp], &z[0], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], info);
//Apply transformation matrix used in reduction to tridiagonal
//form to eigenvectors returned by DSTEIN.
	for (j = 0; j < (*m); j++) {
	    Rcopy(n, &z[j * ldz + 1], 1, &work[0], 1);
	    Rgemv("N", n, n, One, &q[0], ldq, &work[0], 1, Zero, &z[j * ldz + 1], 1);
	}
    }
  L30:
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
		Rswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
		if (*info != 0) {
		    itmp1 = ifail[i];
		    ifail[i] = ifail[j];
		    ifail[j] = itmp1;
		}
	    }

	}
    }
    return;
}
