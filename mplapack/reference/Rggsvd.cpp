/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggsvd.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rggsvd(const char *jobu, const char *jobv, const char *jobq, INTEGER m,
	    INTEGER n, INTEGER p, INTEGER * k, INTEGER * l, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * alpha, REAL * beta,
	    REAL * u, INTEGER ldu, REAL * v, INTEGER ldv, REAL * q, INTEGER ldq, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j;
    REAL ulp;
    INTEGER ibnd;
    REAL tola;
    INTEGER isub;
    REAL tolb, unfl, temp, smax;
    REAL anorm, bnorm;
    INTEGER wantq, wantu, wantv, ncycle;

    ncycle = 0;
//Test the input parameters
    wantu = Mlsame(jobu, "U");
    wantv = Mlsame(jobv, "V");
    wantq = Mlsame(jobq, "Q");

    *info = 0;
    if (!(wantu || Mlsame(jobu, "N"))) {
	*info = -1;
    } else if (!(wantv || Mlsame(jobv, "N"))) {
	*info = -2;
    } else if (!(wantq || Mlsame(jobq, "N"))) {
	*info = -3;
    } else if (m < 0) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (p < 0) {
	*info = -6;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -10;
    } else if (ldb < max((INTEGER) 1, p)) {
	*info = -12;
    } else if (ldu < 1 || (wantu && ldu < m)) {
	*info = -16;
    } else if (ldv < 1 || (wantv && ldv < p)) {
	*info = -18;
    } else if (ldq < 1 || (wantq && ldq < n)) {
	*info = -20;
    }
    if (*info != 0) {
	Mxerbla("Rggsvd", -(*info));
	return;
    }
//Compute the Frobenius norm of matrices A and B
    anorm = Rlange("1", m, n, &A[0], lda, &work[0]);
    bnorm = Rlange("1", p, n, &B[0], ldb, &work[0]);
//Get machine precision and set up threshold for determining
//the effective numerical rank of the matrices A and B.
    ulp = Rlamch("Precision");
    unfl = Rlamch("Safe Minimum");
    tola = max(m, n) * max(anorm, unfl) * ulp;
    tolb = max(p, n) * max(bnorm, unfl) * ulp;
//Preprocessing
    Rggsvp(jobu, jobv, jobq, m, p, n, &A[0], lda, &B[0], ldb, tola, tolb, k, l, &u[0], ldu, &v[0], ldv, &q[0], ldq, &iwork[1], &work[0], &work[n + 1], info);
//Compute the GSVD of two upper "triangular" matrices
    Rtgsja(jobu, jobv, jobq, m, p, n, (*k), (*l), &A[0], lda, &B[0], ldb, tola, tolb, &alpha[1], &beta[1], &u[0], ldu, &v[0], ldv, &q[0], ldq, &work[0], ncycle, info);
//Sort the singular values and store the pivot indices in IWORK
//Copy ALPHA to WORK, then sort ALPHA in WORK
    Rcopy(n, &alpha[1], 1, &work[0], 1);
    ibnd = min((*l), m - (*k));
    for (i = 0; i < ibnd; i++) {
//Scan for largest ALPHA(K+I)
	isub = i;
	smax = work[(*k) + i];
	for (j = i + 1; j <= ibnd; j++) {
	    temp = work[(*k) + j];
	    if (temp > smax) {
		isub = j;
		smax = temp;
	    }
	}
	if (isub != i) {
	    work[(*k) + isub] = work[(*k) + i];
	    work[(*k) + i] = smax;
	    iwork[(*k) + i] = (*k) + isub;
	} else {
	    iwork[(*k) + i] = (*k) + i;
	}
    }
    return;
}
