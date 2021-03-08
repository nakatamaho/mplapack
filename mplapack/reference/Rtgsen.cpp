/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgsen.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rtgsen(INTEGER ijob, LOGICAL wantq, LOGICAL wantz,
	    LOGICAL * select, INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * alphar, REAL * alphai, REAL *
	    beta, REAL * q, INTEGER ldq, REAL * z, INTEGER ldz,
	    INTEGER * m, REAL * pl, REAL * pr, REAL * dif, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, k, n1, n2, kk, ks, mn2, ijb;
    REAL eps;
    INTEGER kase;
    LOGICAL pair;
    INTEGER ierr;
    REAL dsum;
    LOGICAL swap;
    INTEGER isave[3];
    LOGICAL wantd;
    INTEGER lwmin, liwmin;
    LOGICAL wantp;
    LOGICAL wantd1, wantd2;
    REAL dscale, rdscal;
    REAL smlnum;
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
	*info = -14;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -16;
    }
    if (*info != 0) {
	Mxerbla("Rtgsen", -(*info));
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    ierr = 0;
    wantp = ijob == 1 || ijob >= 4;
    wantd1 = ijob == 2 || ijob == 4;
    wantd2 = ijob == 3 || ijob == 5;
    wantd = wantd1 || wantd2;
//Set M to the dimension of the specified pair of deflating
//subspaces.
    *m = 0;
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
    if (ijob == 1 || ijob == 2 || ijob == 4) {
	lwmin = max(max((INTEGER) 1, (n * 4) + 16), ((*m) * 2) * (n - (*m)));
	liwmin = max((INTEGER) 1, n + 6);
    } else if (ijob == 3 || ijob == 5) {
	lwmin = max(max((INTEGER) 1, (n << 2) + 16), ((*m) << 2) * (n - (*m)));
	liwmin = max(max((INTEGER) 1, ((*m) << 1) * (n - (*m))), n + 6);
    } else {
	lwmin = max((INTEGER) 1, (n << 2) + 16);
	liwmin = 1;
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    if (lwork < lwmin && !lquery) {
	*info = -22;
    } else if (liwork < liwmin && !lquery) {
	*info = -24;
    }
    if (*info != 0) {
	Mxerbla("Rtgsen", -(*info));
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
		Rlassq(n, &A[i * lda + 1], 1, &dscale, &dsum);
		Rlassq(n, &B[i * ldb + 1], 1, &dscale, &dsum);
	    }
	    dif[1] = dscale * sqrt(dsum);
	    dif[2] = dif[1];
	}
	goto L60;
    }
//Collect the selected blocks at the top-left corner of (A, B).
    ks = 0;
    pair = MFALSE;
    for (k = 0; k < n; k++) {
	if (pair) {
	    pair = MFALSE;
	} else {
	    swap = select[k];
	    if (k < n) {
		if (A[k + 1 + k * lda] != Zero) {
		    pair = MTRUE;
		    swap = swap || select[k + 1];
		}
	    }
	    if (swap) {
		ks++;
//Swap the K-th block to position KS.
//Perform the reordering of diagonal blocks in (A, B)
//by orthogonal transformation matrices and update
//Q and Z accordingly (if requested):
		kk = k;
		if (k != ks) {
		    Rtgexc(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, &kk, &ks, &work[0], lwork, &ierr);
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
		    goto L60;
		}

		if (pair) {
		    ks++;
		}
	    }
	}
    }
    if (wantp) {
//Solve generalized Sylvester equation for R and L
//and compute PL and PR.
	n1 = (*m);
	n2 = n - (*m);
	i = n1 + 1;
	ijb = 0;
	Rlacpy("Full", n1, n2, &A[i * lda + 1], lda, &work[0], n1);
	Rlacpy("Full", n1, n2, &B[i * ldb + 1], ldb, &work[n1 * n2 + 1], n1);
	Rtgsyl("N", ijb, n1, n2, &A[0], lda, &A[i + i * lda], lda, &work[0], n1, &B[0], ldb, &B[i + i * ldb], ldb,
	       &work[n1 * n2 + 1], n1, &dscale, &dif[1], &work[(n1 * n2 << 1) + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
//Estimate the reciprocal of norms of "projections" onto left
//and right eigenspaces.
	rdscal = Zero;
	dsum = One;
	Rlassq(n1 * n2, &work[0], 1, &rdscal, &dsum);
	*pl = rdscal * sqrt(dsum);
	if (*pl == Zero) {
	    *pl = Zero;
	} else {
	    *pl = dscale / (sqrt(dscale * dscale / *pl + *pl) * sqrt(*pl));
	}
	rdscal = Zero;
	dsum = One;
	Rlassq(n1 * n2, &work[n1 * n2 + 1], 1, &rdscal, &dsum);
	*pr = rdscal * sqrt(dsum);
	if (*pr == Zero) {
	    *pr = One;
	} else {
	    *pr = dscale / (sqrt(dscale * dscale / *pr + *pr) * sqrt(*pr));
	}
    }
    if (wantd) {
//Compute estimates of Difu and Difl.
	if (wantd1) {
	    n1 = (*m);
	    n2 = n - (*m);
	    i = n1 + 1;
	    ijb = 3;
//Frobenius norm-based Difu-estimate.
	    Rtgsyl("N", ijb, n1, n2, &A[0], lda, &A[i + i *
						    lda], lda, &work[0], n1, &B[0], ldb, &B[i + i * ldb], ldb, &work[n1 * n2 + 1],
		   n1, &dscale, &dif[1], &work[(n1 << 1) * n2 + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
//Frobenius norm-based Difl-estimate.
	    Rtgsyl("N", ijb, n2, n1, &A[i + i * lda], lda, &A[0], lda, &work[0], n2, &B[i + i * ldb],
		   ldb, &B[0], ldb, &work[n1 * n2 + 1], n2, &dscale, &dif[2], &work[(n1 << 1) * n2 + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
	} else {
//Compute 1-norm-based estimates of Difu and Difl using
//reversed communication with DLACNTwo In each step a
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
	    Rlacn2(mn2, &work[mn2 + 1], &work[0], &iwork[1], &dif[1], &kase, isave);
	    if (kase != 0) {
		if (kase == 1) {
//Solve generalized Sylvester equation.
		    Rtgsyl("N", ijb, n1, n2, &A[0], lda, &A[i + i * lda], lda, &work[0], n1, &B[0],
			   ldb, &B[i + i * ldb], ldb, &work[n1 * n2 + 1], n1, &dscale, &dif[1], &work[(n1 << 1) * n2 + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		} else {
//Solve the transposed variant.
		    Rtgsyl("T", ijb, n1, n2, &A[0], lda, &A[i + i * lda], lda, &work[0], n1, &B[0],
			   ldb, &B[i + i * ldb], ldb, &work[n1 * n2 + 1], n1, &dscale, &dif[1], &work[(n1 << 1) * n2 + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		}
		goto L40;
	    }
	    dif[1] = dscale / dif[1];
//1-norm-based estimate of Difl.
	  L50:
	    Rlacn2(mn2, &work[mn2 + 1], &work[0], &iwork[1], &dif[2], &kase, isave);
	    if (kase != 0) {
		if (kase == 1) {
//Solve generalized Sylvester equation.
		    Rtgsyl("N", ijb, n2, n1, &A[i + i * lda], lda,
			   &A[0], lda, &work[0], n2, &B[i + i * ldb], ldb, &B[0], ldb, &work[n1 * n2 + 1], n2, &dscale, &dif[2],
			   &work[(n1 << 1) * n2 + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		} else {
//Solve the transposed variant.
		    Rtgsyl("T", ijb, n2, n1, &A[i + i * lda], lda,
			   &A[0], lda, &work[0], n2, &B[i + i * ldb], ldb, &B[0], ldb, &work[n1 * n2 + 1], n2, &dscale, &dif[2],
			   &work[(n1 << 1) * n2 + 1], lwork - (n1 << 1) * n2, &iwork[1], &ierr);
		}
		goto L50;
	    }
	    dif[2] = dscale / dif[2];
	}
    }
  L60:
//Compute generalized eigenvalues of reordered pair (A, B) and
//normalize the generalized Schur form.
    pair = MFALSE;
    for (k = 0; k < n; k++) {
	if (pair) {
	    pair = MFALSE;
	} else {
	    if (k < n) {
		if (A[k + 1 + k * lda] != Zero) {
		    pair = MTRUE;
		}
	    }
	    if (pair) {
//Compute the eigenvalue(s) at position K.
		work[1] = A[k + k * lda];
		work[2] = A[k + 1 + k * lda];
		work[3] = A[k + (k + 1) * lda];
		work[4] = A[k + 1 + (k + 1) * lda];
		work[5] = B[k + k * ldb];
		work[6] = B[k + 1 + k * ldb];
		work[7] = B[k + (k + 1) * ldb];
		work[8] = B[k + 1 + (k + 1) * ldb];
		Rlag2(&work[0], 2, &work[5], 2, smlnum * eps, &beta[k], &beta[k + 1], &alphar[k], &alphar[k + 1], &alphai[k]);
		alphai[k + 1] = -alphai[k];
	    } else {
		if (sign(One, B[k + k * ldb]) < Zero) {
//If B(K,K) is negative, make it positive
		    for (i = 0; i < n; i++) {
			A[k + i * lda] = -A[k + i * lda];
			B[k + i * ldb] = -B[k + i * ldb];
			q[i + k * ldq] = -q[i + k * ldq];
		    }
		}
		alphar[k] = A[k + k * lda];
		alphai[k] = Zero;
		beta[k] = B[k + k * ldb];
	    }
	}
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
