/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlagv2.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlagv2(REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * alphar, REAL * alphai, REAL * beta, REAL * csl, REAL * snl, REAL * csr, REAL * snr)
{
    REAL r, t, h1, h2, h3, wi, qq, rr, wr1, wr2, ulp;
    REAL ascale, bscale, anorm, bnorm;
    REAL safmin, scale1, scale2;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2, mtemp3;

    safmin = Rlamch("S");
    ulp = Rlamch("P");

//Scale A
    mtemp1 = abs(A[lda + 1]) + abs(A[lda + 2]), mtemp2 = abs(A[(lda * 2) + 1]) + abs(A[(lda * 2) + 2]);
    mtemp3 = max(mtemp1, mtemp2);
    anorm = max(mtemp3, safmin);
    ascale = One / anorm;
    A[lda + 1] = ascale * A[lda + 1];
    A[lda * 2 + 1] = ascale * A[lda * 2 + 1];
    A[lda + 2] = ascale * A[lda + 2];
    A[lda * 2 + 2] = ascale * A[lda * 2 + 2];

//Scale B
    mtemp1 = abs(B[ldb + 1]), mtemp2 = abs(B[(ldb * 2) + 1]) + abs(B[(ldb * 2) + 2]);
    mtemp3 = max(mtemp1, mtemp2);
    bnorm = max(mtemp3, safmin);
    bscale = One / bnorm;
    B[ldb + 1] = bscale * B[ldb + 1];
    B[(ldb * 2) + 1] = bscale * B[(ldb * 2) + 1];
    B[(ldb * 2) + 2] = bscale * B[(ldb * 2) + 2];
//Check if A can be deflated
    if (abs(A[lda + 2]) <= ulp) {
	*csl = Zero;
	*snl = Zero;
	*csr = One;
	*snr = Zero;
	A[lda + 2] = Zero;
	B[ldb + 2] = Zero;
//Check if B is singular
    } else if (abs(B[ldb + 1]) <= ulp) {
	Rlartg(A[lda + 1], A[lda + 2], csl, snl, &r);
	*csr = One;
	*snr = Zero;
	Rrot(2, &A[lda + 1], lda, &A[lda + 2], lda, *csl, *snl);
	Rrot(2, &B[ldb + 1], ldb, &B[ldb + 2], ldb, *csl, *snl);
	A[lda + 2] = Zero;
	B[ldb + 1] = Zero;
	B[ldb + 2] = Zero;

    } else if (abs(B[(ldb * 2) + 2]) <= ulp) {
	Rlartg(A[lda * 2 + 2], A[lda + 2], csr, snr, &t);
	*snr = -(*snr);
	Rrot(2, &A[lda + 1], 1, &A[lda * 2 + 1], 1, *csr, *snr);
	Rrot(2, &B[ldb + 1], 1, &B[(ldb * 2) + 1], 1, *csr, *snr);
	*csl = Zero;
	*snl = Zero;
	A[lda + 2] = Zero;
	B[ldb + 2] = Zero;
	B[(ldb * 2) + 2] = Zero;

    } else {
//B is nonsingular, first compute the eigenvalues of (A,B)
	Rlag2(&A[0], lda, &B[0], ldb, safmin, &scale1, &scale2, &wr1, &wr2, &wi);
	if (wi == Zero) {
//two real eigenvalues, compute s*A-w*B
	    h1 = scale1 * A[lda + 1] - wr1 * B[ldb + 1];
	    h2 = scale1 * A[lda * 2 + 1] - wr1 * B[(ldb * 2) + 1];
	    h3 = scale1 * A[lda * 2 + 2] - wr1 * B[(ldb * 2) + 2];
	    rr = Rlapy2(h1, h2);
	    mtemp1 = scale1 * A[lda + 2];
	    qq = Rlapy2(mtemp1, h3);
	    if (rr > qq) {
//find right rotation matrix to zero 1,1 element of
//(sA - wB)
		Rlartg(h2, h1, csr, snr, &t);
	    } else {
//find right rotation matrix to zero 2,1 element of
//(sA - wB)
		Rlartg(h3, scale1 * A[lda + 2], csr, snr, &t);
	    }
	    *snr = -(*snr);
	    Rrot(2, &A[lda + 1], 1, &A[lda * 2 + 1], 1, *csr, *snr);
	    Rrot(2, &B[ldb + 1], 1, &B[(ldb * 2) + 1], 1, *csr, *snr);
//compute inf norms of A and B
	    mtemp1 = abs(A[lda + 1]) + abs(A[lda * 2 + 1]);
	    mtemp2 = abs(A[lda + 2]) + abs(A[lda * 2 + 2]);
	    h1 = max(mtemp1, mtemp2);
	    mtemp1 = abs(B[ldb + 1]) + abs(B[(ldb * 2) + 1]);
	    mtemp2 = abs(B[ldb + 2]) + abs(B[(ldb * 2) + 2]);
	    h2 = max(mtemp1, mtemp2);
	    if (scale1 * h1 >= abs(wr1) * h2) {
//find left rotation matrix Q to zero out B(2,1)
		Rlartg(B[ldb + 1], B[ldb + 2], csl, snl, &r);
	    } else {
//find left rotation matrix Q to zero out A(2,1)
		Rlartg(A[lda + 1], A[lda + 2], csl, snl, &r);
	    }
	    Rrot(2, &A[lda + 1], lda, &A[lda + 2], lda, *csl, *snl);
	    Rrot(2, &B[ldb + 1], ldb, &B[ldb + 2], ldb, *csl, *snl);
	    A[lda + 2] = Zero;
	    B[ldb + 2] = Zero;
	} else {
//a pair of complex conjugate eigenvalues
//first compute the SVD of the matrix B
	    Rlasv2(B[ldb + 1], B[(ldb * 2) + 1], B[(ldb * 2) + 2], &r, &t, snr, csr, snl, csl);
//Form (A,B) := Q(A,B)Z' where Q is left rotation matrix and
//Z is right rotation matrix computed from DLASV2

	    Rrot(2, &A[lda + 1], lda, &A[lda + 2], lda, *csl, *snl);
	    Rrot(2, &B[ldb + 1], ldb, &B[ldb + 2], ldb, *csl, *snl);
	    Rrot(2, &A[lda + 1], 1, &A[lda * 2 + 1], 1, *csr, *snr);
	    Rrot(2, &B[ldb + 1], 1, &B[(ldb * 2) + 1], 1, *csr, *snr);
	    B[ldb + 2] = Zero;
	    B[(ldb * 2) + 1] = Zero;
	}
    }
//Unscaling
    A[lda + 1] = anorm * A[lda + 1];
    A[lda + 2] = anorm * A[lda + 2];
    A[lda * 2 + 1] = anorm * A[lda * 2 + 1];
    A[lda * 2 + 2] = anorm * A[lda * 2 + 2];
    B[ldb + 1] = bnorm * B[ldb + 1];
    B[ldb + 2] = bnorm * B[ldb + 2];
    B[(ldb * 2) + 1] = bnorm * B[(ldb * 2) + 1];
    B[(ldb * 2) + 2] = bnorm * B[(ldb * 2) + 2];

    if (wi == Zero) {
	alphar[1] = A[lda + 1];
	alphar[2] = A[lda * 2 + 2];
	alphai[1] = Zero;
	alphai[2] = Zero;
	beta[1] = B[ldb + 1];
	beta[2] = B[(ldb * 2) + 2];
    } else {
	alphar[1] = anorm * wr1 / scale1 / bnorm;
	alphai[1] = anorm * wi / scale1 / bnorm;
	alphar[2] = alphar[1];
	alphai[2] = -alphai[1];
	beta[1] = One;
	beta[2] = One;
    }
    return;
}
