/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgsna.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rtgsna(const char *job, const char *howmny, LOGICAL * select,
	    INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb,
	    REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, REAL * s, REAL * dif, INTEGER mm, INTEGER * m, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, k;
    REAL c1, c2;
    INTEGER n1, n2, ks, iz;
    REAL eps, beta, cond = 0;
    LOGICAL pair;
    INTEGER ierr;
    REAL uhav, uhbv;
    INTEGER ifst;
    REAL lnrm;
    INTEGER ilst;
    REAL rnrm;
    REAL root1, root2, scale;
    REAL uhavi, uhbvi, tmpii;
    INTEGER lwmin;
    LOGICAL wants;
    REAL tmpir, tmpri, dummy[1], tmprr;
    REAL dummy1[1];
    REAL alphai, alphar;
    LOGICAL wantbh, wantdf, somcon;
    REAL alprqt = 0.0;
    REAL smlnum;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Four = 4.0;
    REAL mtemp1, mtemp2;

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
	    (*m) = 0;
	    pair = MFALSE;
	    for (k = 0; k < n; k++) {
		if (pair) {
		    pair = MFALSE;
		} else {
		    if (k < n) {
			if (A[k + 1 + k * lda] == Zero) {
			    if (select[k]) {
				++(*m);
			    }
			} else {
			    pair = MTRUE;
			    if (select[k] || select[k + 1]) {
				(*m) = (*m) + 2;
			    }
			}
		    } else {
			if (select[n]) {
			    ++(*m);
			}
		    }
		}
	    }
	} else {
	    (*m) = n;
	}

	if (n == 0) {
	    lwmin = 1;
	} else if (Mlsame(job, "V") || Mlsame(job, "B")) {
	    lwmin = (n << 1) * (n + 2) + 16;
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
	Mxerbla("Rtgsna", -(*info));
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
    ks = 0;
    pair = MFALSE;
    for (k = 0; k < n; k++) {
//Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block.
	if (pair) {
	    pair = MFALSE;
	    goto L20;
	} else {
	    if (k < n) {
		pair = A[k + 1 + k * lda] != Zero;
	    }
	}
//Determine whether condition numbers are required for the k-th
//eigenpair.
	if (somcon) {
	    if (pair) {
		if (!select[k] && !select[k + 1]) {
		    goto L20;
		}
	    } else {
		if (!select[k]) {
		    goto L20;
		}
	    }
	}
	ks++;
	if (wants) {
//Compute the reciprocal condition number of the k-th
//eigenvalue.
	    if (pair) {
//Complex eigenvalue pair.
		mtemp1 = Rnrm2(n, &vr[ks * ldvr + 1], 1);
		mtemp2 = Rnrm2(n, &vr[(ks + 1) * ldvr + 1], 1);
		rnrm = Rlapy2(mtemp1, mtemp2);
		mtemp1 = Rnrm2(n, &vl[ks * ldvl + 1], 1);
		mtemp2 = Rnrm2(n, &vl[(ks + 1) * ldvl + 1], 1);
		lnrm = Rlapy2(mtemp1, mtemp2);
		Rgemv("N", n, n, One, &A[0], lda, &vr[ks * ldvr + 1], 1, Zero, &work[0], 1);
		tmprr = Rdot(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
		tmpri = Rdot(n, &work[0], 1, &vl[(ks + 1) * ldvl + 1], 1);
		Rgemv("N", n, n, One, &A[0], lda, &vr[(ks + 1) * ldvr + 1], 1, Zero, &work[0], 1);
		tmpii = Rdot(n, &work[0], 1, &vl[(ks + 1) * ldvl + 1], 1);
		tmpir = Rdot(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
		uhav = tmprr + tmpii;
		uhavi = tmpir - tmpri;
		Rgemv("N", n, n, One, &B[0], ldb, &vr[ks * ldvr + 1], 1, Zero, &work[0], 1);
		tmprr = Rdot(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
		tmpri = Rdot(n, &work[0], 1, &vl[(ks + 1) * ldvl + 1], 1);
		Rgemv("N", n, n, One, &B[0], ldb, &vr[(ks + 1) * ldvr + 1], 1, Zero, &work[0], 1);
		tmpii = Rdot(n, &work[0], 1, &vl[(ks + 1) * ldvl + 1], 1);
		tmpir = Rdot(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
		uhbv = tmprr + tmpii;
		uhbvi = tmpir - tmpri;
		uhav = Rlapy2(uhav, uhavi);
		uhbv = Rlapy2(uhbv, uhbvi);
		cond = Rlapy2(uhav, uhbv);
		s[ks] = cond / (rnrm * lnrm);
		s[ks + 1] = s[ks];

	    } else {
//Real eigenvalue.
		rnrm = Rnrm2(n, &vr[ks * ldvr + 1], 1);
		lnrm = Rnrm2(n, &vl[ks * ldvl + 1], 1);
		Rgemv("N", n, n, One, &A[0], lda, &vr[ks * ldvr + 1], 1, Zero, &work[0], 1);
		uhav = Rdot(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
		Rgemv("N", n, n, One, &B[0], ldb, &vr[ks * ldvr + 1], 1, Zero, &work[0], 1);
		uhbv = Rdot(n, &work[0], 1, &vl[ks * ldvl + 1], 1);
		cond = Rlapy2(uhav, uhbv);
		if (cond == Zero) {
		    s[ks] = -One;
		} else {
		    s[ks] = cond / (rnrm * lnrm);
		}
	    }
	}
	if (wantdf) {
	    if (n == 1) {
		dif[ks] = Rlapy2(A[lda + 1], B[ldb + 1]);
		goto L20;
	    }
//Estimate the reciprocal condition number of the k-th
//eigenvectors.
	    if (pair) {
//Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)).
//Compute the eigenvalue(s) at position K.
		work[1] = A[k + k * lda];
		work[2] = A[k + 1 + k * lda];
		work[3] = A[k + (k + 1) * lda];
		work[4] = A[k + 1 + (k + 1) * lda];
		work[5] = B[k + k * ldb];
		work[6] = B[k + 1 + k * ldb];
		work[7] = B[k + (k + 1) * ldb];
		work[8] = B[k + 1 + (k + 1) * ldb];
		Rlag2(&work[0], 2, &work[5], 2, smlnum * eps, &beta, dummy1, &alphar, dummy, &alphai);
		alprqt = One;
		c1 = (alphar * alphar + alphai * alphai + beta * beta) * Two;
		c2 = beta * Four * beta * alphai * alphai;
		root1 = c1 + sqrt(c1 * c1 - c2 * Four);
		root2 = c2 / root1;
		root1 = root1 / Two;
		mtemp1 = sqrt(root1), mtemp2 = sqrt(root2);
		cond = min(mtemp1, mtemp2);
	    }
//Copy the matrix (A, B) to the array WORK and swap the
//diagonal block beginning at A(k,k) to the (1,1) position.
	    Rlacpy("Full", n, n, &A[0], lda, &work[0], n);
	    Rlacpy("Full", n, n, &B[0], ldb, &work[n * n + 1], n);
	    ifst = k;
	    ilst = 1;
	    Rtgexc(MFALSE, MFALSE, n, &work[0], n, &work[n * n + 1], n, dummy, 1, dummy1, 1, &ifst, &ilst, &work[(n * n << 1) + 1], lwork - (n << 1) * n, &ierr);
	    if (ierr > 0) {
//Ill-conditioned problem - swap rejected.
		dif[ks] = Zero;
	    } else {
//Reordering successful, solve generalized Sylvester
//equation for R and L,
//           A22 * R - L * A11 = A12
//           B22 * R - L * B11 = B12,
//and compute estimate of Difl((A11,B11), (A22, B22)).
		n1 = 1;
		if (work[2] != Zero) {
		    n1 = 2;
		}
		n2 = n - n1;
		if (n2 == 0) {
		    dif[ks] = cond;
		} else {
		    i = n * n + 1;
		    iz = (n << 1) * n + 1;
		    Rtgsyl("N", 3, n2, n1, &work[n * n1 + n1 + 1], n,
			   &work[0], n, &work[n1 + 1], n, &work[n * n1 + n1 + i], n, &work[i], n, &work[n1 + i], n, &scale,
			   &dif[ks], &work[iz + 1], lwork - (n << 1) * n, &iwork[1], &ierr);
		    if (pair) {
			mtemp1 = max(One, alprqt) * dif[ks];
			dif[ks] = min(mtemp1, cond);
		    }
		}
	    }
	    if (pair) {
		dif[ks + 1] = dif[ks];
	    }
	}
	if (pair) {
	    ks++;
	}
      L20:
	;
    }
    work[1] = lwmin;
    return;
}
