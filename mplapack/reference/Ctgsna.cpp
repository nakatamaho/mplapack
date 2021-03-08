/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctgsna.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MFALSE 0
#define MTRUE 1

void Ctgsna(const char *job, const char *howmny, LOGICAL * select, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb,
	    COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, REAL * s, REAL * dif, INTEGER mm, INTEGER * m, COMPLEX * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, k, n1 = 1, ks;
    REAL eps, cond;
    INTEGER ierr, ifst;
    REAL lnrm;
    COMPLEX yhax, yhbx;
    INTEGER ilst;
    REAL rnrm, scale;
    INTEGER lwmin;
    INTEGER wants;
    COMPLEX dummy[1];
    COMPLEX dummy1[1];
    REAL bignum;
    INTEGER wantbh, wantdf, somcon;
    REAL smlnum;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;

//Decode and test the input parameters
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantdf = Mlsame(job, "V") || wantbh;
    somcon = Mlsame(howmny, "S");
    *info = 0;
    lquery = lwork == -1;
    if (!wants && !wantdf) {
	*info = -1;
    } else if (!Mlsame(howmny, "A") && !somcon) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    } else if (wants && ldvl < n) {
	*info = -10;
    } else if (wants && ldvr < n) {
	*info = -12;
    } else {
//Set M to the number of eigenpairs for which condition numbers
//are required, and test MM.
	if (somcon) {
	    m = 0;
	    for (k = 0; k < n; k++) {
		if (select[k]) {
		    ++(*m);
		}
	    }
	} else {
	    (*m) = n;
	}
	if (n == 0) {
	    lwmin = 1;
	} else if (Mlsame(job, "V") || Mlsame(job, "B")) {
	    lwmin = (n << 1) * n;
	} else {
	    lwmin = n;
	}
	work[1] = lwmin;
	if (mm < (*m)) {
	    *info = -15;
	} else if (lwork < lwmin && !lquery) {
	    *info = -18;
	}
    }
    if (*info != 0) {
	Mxerbla("Ctgsna", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = One / smlnum;
    ks = 0;
    for (k = 0; k < n; k++) {
//Determine whether condition numbers are required for the k-th
//eigenpair.
	if (somcon) {
	    if (!select[k]) {
		goto L20;
	    }
	}
	ks++;
	if (wants) {
//Compute the reciprocal condition number of the k-th
//eigenvalue.
	    rnrm = RCnrm2(n, &vr[ks * ldvr + 1], 1);
	    lnrm = RCnrm2(n, &vl[ks * ldvl + 1], 1);
	    Cgemv("N", n, n, (COMPLEX) One, &A[0], lda, &vr[ks * ldvr + 1], 1, (COMPLEX) Zero, &work[0], 1);
	    yhax = Cdotc(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
	    Cgemv("N", n, n, (COMPLEX) One, &B[0], ldb, &vr[ks * ldvr + 1], 1, (COMPLEX) Zero, &work[0], 1);
	    yhbx = Cdotc(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
	    cond = Rlapy2(abs(yhax), abs(yhbx));
	    if (cond == Zero) {
		s[ks] = -One;
	    } else {
		s[ks] = cond / (rnrm * lnrm);
	    }
	}
	if (wantdf) {
	    if (n == 1) {
		dif[ks] = Rlapy2(abs(A[lda + 1]), abs(B[ldb + 1]));
	    } else {
//Estimate the reciprocal condition number of the k-th
//eigenvectors.
//Copy the matrix (A, B) to the array WORK and move the
//(k,k)th pair to the (1,1) position.
		Clacpy("Full", n, n, &A[0], lda, &work[0], n);
		Clacpy("Full", n, n, &B[0], ldb, &work[n * n + 1], n);
		ifst = k;
		ilst = 1;
		Ctgexc(MFALSE, MFALSE, n, &work[0], n, &work[n * n + 1], n, dummy, 1, dummy1, 1, ifst, &ilst, &ierr);
		if (ierr > 0) {
//Ill-conditioned problem - swap rejected.
		    dif[ks] = Zero;
		} else {
//Reordering successful, solve generalized Sylvester */
//equation for R and L, */
//           A22 * R - L * A11 = A12 */
//           B22 * R - L * B11 = B12, */
//and compute estimate of Difl[(A11,B11), (A22, B22)]. */
		    i = n * n + 1;
		    Ctgsyl("N", 3, n - n1, 1, &work[n * n1 + n1 + 1], n,
			   &work[0], n, &work[n1 + 1], n, &work[n * n1 + n1 + i], n, &work[i], n, &work[n1 + i], n, &scale, &dif[ks], dummy, 1, &iwork[1], &ierr);
		}
	    }
	}
      L20:
	;
    }
    work[1] = lwmin;
    return;
}
