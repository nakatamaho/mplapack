/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cggsvp.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cggsvp(const char *jobu, const char *jobv, const char *jobq, INTEGER m,
	    INTEGER p, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, REAL tola, REAL tolb, INTEGER * k,
	    INTEGER * l, COMPLEX * u, INTEGER ldu, COMPLEX * v, INTEGER ldv, COMPLEX * q, INTEGER ldq, INTEGER * iwork, REAL * rwork, COMPLEX * tau, COMPLEX * work, INTEGER * info)
{
    INTEGER i, j;
    INTEGER wantq, wantu, wantv;
    INTEGER forwrd;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters
    wantu = Mlsame(jobu, "U");
    wantv = Mlsame(jobv, "V");
    wantq = Mlsame(jobq, "Q");
    forwrd = MTRUE;
    *info = 0;
    if (!(wantu || Mlsame(jobu, "N"))) {
	*info = -1;
    } else if (!(wantv || Mlsame(jobv, "N"))) {
	*info = -2;
    } else if (!(wantq || Mlsame(jobq, "N"))) {
	*info = -3;
    } else if (m < 0) {
	*info = -4;
    } else if (p < 0) {
	*info = -5;
    } else if (n < 0) {
	*info = -6;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -8;
    } else if (ldb < max((INTEGER) 1, p)) {
	*info = -10;
    } else if (ldu < 1 || (wantu && ldu < m)) {
	*info = -16;
    } else if (ldv < 1 || (wantv && ldv < p)) {
	*info = -18;
    } else if (ldq < 1 || (wantq && ldq < n)) {
	*info = -20;
    }
    if (*info != 0) {
	Mxerbla("Cggsvp", -(*info));
	return;
    }
//QR with column pivoting of B: B*P = V*( S11 S12 )
//                                      (  0   0  )
    for (i = 0; i < n; i++) {
	iwork[i] = 0;
    }
    Cgeqpf(p, n, &B[0], ldb, &iwork[1], &tau[1], &work[0], &rwork[1], info);
//Update A := A*P
    Clapmt(&forwrd, m, n, &A[0], lda, &iwork[1]);
//Determine the effective rank of matrix B.
    l = 0;
    for (i = 0; i < min(p, n); i++) {
	if (Cabs1(B[i + i * ldb]) > tolb) {
	    ++(l);
	}
    }
    if (wantv) {
//Copy the details of V, and form V.
	Claset("Full", p, p, Zero, Zero, &v[0], ldv);
	if (p > 1) {
	    Clacpy("Lower", p - 1, n, &B[ldb + 2], ldb, &v[ldv + 2], ldv);
	}
	Cung2r(p, p, min(p, n), &v[0], ldv, &tau[1], &work[0], info);
    }
//Clean up B
    for (j = 0; j < (*l) - 1; j++) {
	for (i = j + 1; i <= (*l); i++) {
	    B[i + j * ldb] = Zero;
	}
    }
    if (p > (*l)) {
	Claset("Full", p - (*l), n, Zero, Zero, &B[(*l) + 1 + ldb], ldb);
    }
    if (wantq) {
//Set Q = I and Update Q := Q*P
	Claset("Full", n, n, Zero, One, &q[0], ldq);
	Clapmt(&forwrd, n, n, &q[0], ldq, &iwork[1]);
    }
    if (p >= (*l) && n != (*l)) {
//RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z
	Cgerq2((*l), n, &B[0], ldb, &tau[1], &work[0], info);
//Update A := A*Z'
	Cunmr2("Right", "Conjugate transpose", m, n, (*l), &B[0], ldb, &tau[1], &A[0], lda, &work[0], info);
	if (wantq) {
//Update Q := Q*Z'
	    Cunmr2("Right", "Conjugate transpose", n, n, (*l), &B[0], ldb, &tau[1], &q[0], ldq, &work[0], info);
	}
//Clean up B
	Claset("Full", (*l), n - (*l), Zero, Zero, &B[0], ldb);
	for (j = n - (*l) + 1; j <= n; j++) {
	    for (i = j - n + (*l) + 1; i <= (*l); i++) {
		B[i + j * ldb] = Zero;
	    }
	}
    }
//Let              N-L     L
//           A = ( A11    A12 ) M,
//then the following does the complete QR decomposition of A11:
//         A11 = U*(  0  T12 )*P1'
//                 (  0   0  )
    for (i = 0; i < n - (*l); i++) {
	iwork[i] = 0;
    }
    Cgeqpf(m, n - (*l), &A[0], lda, &iwork[1], &tau[1], &work[0], &rwork[1], info);
//Determine the effective rank of A11
    k = 0;
    for (i = 0; i < min(m, n - (*l)); i++) {
	if (Cabs1(A[i + i * lda]) > tola) {
	    ++(k);
	}
    }
//Update A12 := U'*A12, where A12 = A( 1:M, N-L+1:N )
    Cunm2r("Left", "Conjugate transpose", m, (*l), min(m, n - (*l)), &A[0], lda, &tau[1], &A[(n - (*l) + 1) * lda], lda, &work[0], info);
    if (wantu) {
//Copy the details of U, and form U
	Claset("Full", m, m, Zero, Zero, &u[0], ldu);
	if (m > 1) {
	    Clacpy("Lower", m - 1, n - (*l), &A[lda + 2], lda, &u[ldu + 2]
		   , ldu);
	}
	Cung2r(m, m, min(m, n - (*l)), &u[0], ldu, &tau[1], &work[0], info);
    }
    if (wantq) {
//Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
	Clapmt(&forwrd, n, n - (*l), &q[0], ldq, &iwork[1]);
    }
//Clean up A: set the strictly lower triangular part of
//A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = Zero
    for (j = 0; j < (*k) - 1; j++) {
	for (i = j + 1; i <= (*k); i++) {
	    A[i + j * lda] = Zero;
	}
    }
    if (m > (*k)) {
	Claset("Full", m - (*k), n - (*l), Zero, Zero, &A[(*k) + 1 + lda], lda);
    }
    if (n - (*l) > (*k)) {
//RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
	Cgerq2((*k), n - (*l), &A[0], lda, &tau[1], &work[0], info);
	if (wantq) {
//Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1'
	    Cunmr2("Right", "Conjugate transpose", n, n - (*l), (*k), &A[0], lda, &tau[1], &q[0], ldq, &work[0], info);
	}
//Clean up A
	Claset("Full", (*k), n - (*l) - (*k), Zero, Zero, &A[0], lda);
	for (j = n - (*l) - (*k) + 1; j <= n - (*l); j++) {
	    for (i = j - n + (*l) + (*k) + 1; i <= (*k); i++) {
		A[i + j * lda] = Zero;
	    }
	}
    }
    if (m > (*k)) {
//QR factorization of A( K+1:M,N-L+1:N )
	Cgeqr2(m - (*k), (*l), &A[(*k) + 1 + (n - (*l) + 1) * lda], lda, &tau[1], &work[0], info);
	if (wantu) {
//Update U(:,K+1:M) := U(:,K+1:M)*U1
	    Cunm2r("Right", "No transpose", m, m - (*k), min(m - (*k), (*l)), &A[(*k) + 1 + (n - (*l) + 1) * lda], lda, &tau[1], &u[((*k) + 1) * ldu + 1], ldu, &work[0], info);
	}
//Clean up
	for (j = n - (*k) + 1; j <= n; j++) {
	    for (i = j - n + (*k) + (*l) + 1; i <= m; i++) {
		A[i + j * lda] = Zero;
	    }
	}
    }
    return;
}
