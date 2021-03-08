/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctgsen.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctgsen(INTEGER ijob, LOGICAL wantq, LOGICAL wantz,
	    LOGICAL * select, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb,
	    COMPLEX * alpha, COMPLEX * beta, COMPLEX * q, INTEGER ldq, COMPLEX * z, INTEGER ldz,
	    INTEGER * m, REAL * pl, REAL * pr, REAL * dif, COMPLEX * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, k, n1, n2, ks, mn2, ijb, kase, ierr;
    REAL dsum;
    LOGICAL swap;
    COMPLEX temp1, temp2;
    INTEGER isave[3];
    LOGICAL wantd;
    INTEGER lwmin;
    LOGICAL wantp;
    LOGICAL wantd1, wantd2;
    REAL dscale, rdscal, safmin;
    INTEGER liwmin;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0;

//Decode and test the input parameters
    *info = 0;
    lquery = lwork == -1 || liwork == -1;

    if (ijob < 0 || ijob > 5) {
	*info = -1;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else if (ldq < 1 || (wantq && ldq < n)) {
	*info = -13;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -15;
    }
    if (*info != 0) {
	Mxerbla("Ctgsen", -(*info));
	return;
    }
    ierr = 0;
    wantp = ijob == 1 || ijob >= 4;
    wantd1 = ijob == 2 || ijob == 4;
    wantd2 = ijob == 3 || ijob == 5;
    wantd = wantd1 || wantd2;
//Set M to the dimension of the specified pair of deflating
//subspaces.
    (*m) = 0;
    for (k = 0; k < n; k++) {
	alpha[k] = A[k + k * lda];
	beta[k] = B[k + k * ldb];
	if (k < n) {
	    if (select[k]) {
		++(*m);
	    }
	} else {
	    if (select[n]) {
		++(*m);
	    }
	}
    }
    if (ijob == 1 || ijob == 2 || ijob == 4) {
	lwmin = max((INTEGER) 1, ((*m) << 1) * (n - (*m)));
	liwmin = max((INTEGER) 1, n + 2);
    } else if (ijob == 3 || ijob == 5) {
	lwmin = max((INTEGER) 1, ((*m) << 2) * (n - (*m)));
	liwmin = max(max((INTEGER) 1, ((*m) << 1) * (n - (*m))), n + 2);
    } else {
	lwmin = 1;
	liwmin = 1;
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    if (lwork < lwmin && !lquery) {
	*info = -21;
    } else if (liwork < liwmin && !lquery) {
	*info = -23;
    }
    if (*info != 0) {
	Mxerbla("Ctgsen", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible.
    if ((*m) == n || (*m) == 0) {
	if (wantp) {
	    *pl = Zero;
	    *pr = One;
	}
	if (wantd) {
	    dscale = Zero;
	    dsum = One;
	    for (i = 0; i < n; i++) {
		Classq(n, &A[i * lda + 1], 1, &dscale, &dsum);
		Classq(n, &B[i * ldb + 1], 1, &dscale, &dsum);
	    }
	    dif[1] = dscale * sqrt(dsum);
	    dif[2] = dif[1];
	}
	goto L70;
    }
//Get machine constant
    safmin = Rlamch("S");
//Collect the selected blocks at the top-left corner of (A, B).
    ks = 0;
    for (k = 0; k < n; k++) {
	swap = select[k];
	if (swap) {
	    ks++;
//Swap the K-th block to position KS. Compute unitary Q
//and Z that will swap adjacent diagonal blocks in (A, B).
	    if (k != ks) {
		Ctgexc(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, k, &ks, &ierr);
	    }
	    if (ierr > 0) {
//Swap is rejected: exit.
		*info = 1;
		if (wantp) {
		    *pl = Zero;
		    *pr = Zero;
		}
		if (wantd) {
		    dif[1] = Zero;
		    dif[2] = Zero;
		}
		goto L70;
	    }
	}
    }
    if (wantp) {
//Solve generalized Sylvester equation for R and L:
//           A11 * R - L * A22 = A12
//           B11 * R - L * B22 = B12
	n1 = (*m);
	n2 = n - (*m);
	i = n1 + 1;
	Clacpy("Full", n1, n2, &A[i * lda + 1], lda, &work[0], n1);
	Clacpy("Full", n1, n2, &B[i * ldb + 1], ldb, &work[n1 * n2 + 1], n1);
	ijb = 0;
	Ctgsyl("N", ijb, n1, n2, &A[0], lda, &A[i + i * lda]
	       , lda, &work[0], n1, &B[0], ldb, &B[i + i * ldb], ldb, &work[n1 * n2 + 1], n1, &dscale, &dif[1],
	       &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
//Estimate the reciprocal of norms of "projections" onto
//left and right eigenspaces
	rdscal = Zero;
	dsum = One;
	Classq(n1 * n2, &work[0], 1, &rdscal, &dsum);
	*pl = rdscal * sqrt(dsum);
	if (*pl == Zero) {
	    *pl = Zero;
	} else {
	    *pl = dscale / (sqrt(dscale * dscale / *pl + *pl) * sqrt(*pl));
	}
	rdscal = Zero;
	dsum = One;
	Classq(n1 * n2, &work[n1 * n2 + 1], 1, &rdscal, &dsum);
	*pr = rdscal * sqrt(dsum);
	if (*pr == Zero) {
	    *pr = One;
	} else {
	    *pr = dscale / (sqrt(dscale * dscale / *pr + *pr) * sqrt(*pr));
	}
    }
    if (wantd) {
//Compute estimates Difu and Difl.
	if (wantd1) {
	    n1 = (*m);
	    n2 = n - (*m);
	    i = n1 + 1;
	    ijb = 3;
//Frobenius norm-based Difu estimate.
	    Ctgsyl("N", ijb, n1, n2, &A[0], lda, &A[i + i * lda], lda, &work[0], n1, &B[0], ldb, &B[i + i * ldb], ldb,
		   &work[n1 * n2 + 1], n1, &dscale, &dif[1], &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
//Frobenius norm-based Difl estimate.
	    Ctgsyl("N", ijb, n2, n1, &A[i + i * lda], lda, &A[0], lda, &work[0], n2, &B[i + i * ldb],
		   ldb, &B[0], ldb, &work[n1 * n2 + 1], n2, &dscale, &dif[2], &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
	} else {
//Compute 1-norm-based estimates of Difu and Difl using
//reversed communication with ZLACNTwo In each step a
//generalized Sylvester equation or a transposed variant
//is solved.
	    kase = 0;
	    n1 = (*m);
	    n2 = n - (*m);
	    i = n1 + 1;
	    ijb = 0;
	    mn2 = (n1 << 1) * n2;
//1-norm-based estimate of Difu.
	  L40:
	    Clacn2(mn2, &work[mn2 + 1], &work[0], &dif[1], &kase, isave);
	    if (kase != 0) {
		if (kase == 1) {
//Solve generalized Sylvester equation
		    Ctgsyl("N", ijb, n1, n2, &A[0], lda, &A[i + i * lda], lda, &work[0], n1, &B[0],
			   ldb, &B[i + i * ldb], ldb, &work[n1 * n2 + 1], n1, &dscale, &dif[1], &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		} else {
//Solve the transposed variant.
		    Ctgsyl("C", ijb, n1, n2, &A[0], lda, &A[i + i * lda], lda, &work[0], n1, &B[0],
			   ldb, &B[i + i * ldb], ldb, &work[n1 * n2 + 1], n1, &dscale, &dif[1], &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		}
		goto L40;
	    }
	    dif[1] = dscale / dif[1];
//1-norm-based estimate of Difl.
	  L50:
	    Clacn2(mn2, &work[mn2 + 1], &work[0], &dif[2], &kase, isave);
	    if (kase != 0) {
		if (kase == 1) {
//Solve generalized Sylvester equation
		    Ctgsyl("N", ijb, n2, n1, &A[i + i * lda], lda,
			   &A[0], lda, &work[0], n2, &B[i + i * ldb], ldb, &B[0], ldb, &work[n1 * n2 + 1], n2, &dscale, &dif[2],
			   &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		} else {
//Solve the transposed variant.
		    Ctgsyl("C", ijb, n2, n1, &A[i + i * lda], lda, &A[0], lda, &work[0], n2, &B[0], ldb, &B[i + i * ldb], ldb,
			   &work[n1 * n2 + 1], n2, &dscale, &dif[2], &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		}
		goto L50;
	    }
	    dif[2] = dscale / dif[2];
	}
    }
//If B(K,K) is complex, make it real and positive (normalization
//of the generalized Schur form) and Store the generalized
//eigenvalues of reordered pair (A, B)
    for (k = 0; k < n; k++) {
	dscale = abs(B[k + k * ldb]);
	if (dscale > safmin) {
	    temp1 = conj(B[k + k * ldb] / dscale);
	    temp2 = B[k + k * ldb] / dscale;
	    B[k + k * ldb] = dscale;
	    Cscal(n - k, temp1, &B[k + (k + 1) * ldb], ldb);
	    Cscal(n - k + 1, temp1, &A[k + k * lda], lda);
	    if (wantq) {
		Cscal(n, temp2, &q[k * ldq + 1], 1);
	    }
	} else {
	    B[k + k * ldb] = Zero;
	}
	alpha[k] = A[k + k * lda];
	beta[k] = B[k + k * ldb];
    }
  L70:
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
