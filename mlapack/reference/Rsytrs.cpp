/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsytrs.cpp,v 1.10 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>

void Rsytrs(const char *uplo, INTEGER n, INTEGER nrhs, REAL * A, INTEGER lda, INTEGER * ipiv, REAL * B, INTEGER ldb, INTEGER * info)
{

    INTEGER j, k;
    REAL ak, bk;
    INTEGER kp;
    REAL akm1, bkm1;
    REAL akm1k;
    REAL denom;
    INTEGER upper;
    REAL One = 1.0, mtemp;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (nrhs < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Rsytrs", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0 || nrhs == 0)
	return;

    if (upper) {
// Solve A*X = B, where A = U*D*U'.
// First solve U*D*X = B, overwriting B with X.
// K is the main loop index, decreasing from N to 1 in steps of
// 1 or 2, depending on the size of the diagonal blocks.

	k = n;
	while (1) {
//If K < 1, exit from loop.
	    if (k < 1)
		break;
	    if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Interchange rows K and IPIV(K).
		kp = ipiv[k];
		if (kp != k) {
		    Rswap(nrhs, &B[k + ldb], ldb, &B[kp + ldb], ldb);
		}
//Multiply by inv(U(K)), where U(K) is the transformation
//stored in column K of A.
		Rger(k - 1, nrhs, -One, &A[k * lda], 1, &B[k + ldb], ldb, &B[ldb + 1], ldb);

//Multiply by the inverse of the diagonal block.
		mtemp = One / A[k + k * lda];
		Rscal(nrhs, mtemp, &B[k + ldb], ldb);
		k--;
	    } else {
//2 x 2 diagonal block
//Interchange rows K-1 and -IPIV(K).
		kp = -ipiv[k];
		if (kp != k - 1) {
		    Rswap(nrhs, &B[k - 1 + ldb], ldb, &B[kp + ldb], ldb);
		}
//Multiply by inv(U(K)), where U(K) is the transformation
//stored in columns K-1 and K of A.
		Rger(k - 2, nrhs, -One, &A[k * lda], 1, &B[k + ldb], ldb, &B[ldb + 1], ldb);
		Rger(k - 2, nrhs, -One, &A[(k - 1) * lda], 1, &B[k - 1 + ldb], ldb, &B[ldb + 1], ldb);

//Multiply by the inverse of the diagonal block.
		akm1k = A[k - 1 + k * lda];
		akm1 = A[k - 1 + (k - 1) * lda] / akm1k;
		ak = A[k + k * lda] / akm1k;
		denom = akm1 * ak - One;
		for (j = 0; j < nrhs; j++) {
		    bkm1 = B[k - 1 + j * ldb] / akm1k;
		    bk = B[k + j * ldb] / akm1k;
		    B[k - 1 + j * ldb] = (ak * bkm1 - bk) / denom;
		    B[k + j * ldb] = (akm1 * bk - bkm1) / denom;

		}
		k += -2;
	    }
	}

//Next solve U'*X = B, overwriting B with X.
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2, depending on the size of the diagonal blocks.

	k = 0;
	while (1) {
//If K > N, exit from loop.
	    if (k > n)
		break;
	    if (ipiv[k] > 0) {

//1 x 1 diagonal block
///Multiply by inv(U'(K)), where U(K) is the transformation
//stored in column K of A.
		Rgemv("Transpose", k - 1, nrhs, -One, &B[0], ldb, &A[k * lda], 1, One, &B[k + ldb], ldb);

//Interchange rows K and IPIV(K).
		kp = ipiv[k];
		if (kp != k) {
		    Rswap(nrhs, &B[k + ldb], ldb, &B[kp + ldb], ldb);
		}
		k++;
	    } else {

//2 x 2 diagonal block
//Multiply by inv(U'(K+1)), where U(K+1) is the transformation
//stored in columns K and K+1 of A.

		Rgemv("Transpose", k - 1, nrhs, -One, &B[0], ldb, &A[k * lda], 1, One, &B[k + ldb], ldb);
		Rgemv("Transpose", k - 1, nrhs, -One, &B[0], ldb, &A[(k + 1) * lda], 1, One, &B[k + 1 + ldb], ldb);

//Interchange rows K and -IPIV(K).
		kp = -ipiv[k];
		if (kp != k) {
		    Rswap(nrhs, &B[k + ldb], ldb, &B[kp + ldb], ldb);
		}
		k += 2;
	    }
	}
    } else {
//Solve A*X = B, where A = L*D*L'.
//First solve L*D*X = B, overwriting B with X.
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2, depending on the size of the diagonal blocks.

	k = 0;
	while (1) {
//If K > N, exit from loop.
	    if (k > n)
		break;
	    if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Interchange rows K and IPIV(K).
		kp = ipiv[k];
		if (kp != k) {
		    Rswap(nrhs, &B[k + ldb], ldb, &B[kp + ldb], ldb);
		}
//Multiply by inv(L(K)), where L(K) is the transformation
//stored in column K of A.
		if (k < n) {
		    Rger(n - k, nrhs, -One, &A[k + 1 + k * lda], 1, &B[k + ldb], ldb, &B[k + 1 + ldb], ldb);
		}
//Multiply by the inverse of the diagonal block.

		mtemp = One / A[k + k * lda];
		Rscal(nrhs, mtemp, &B[k + ldb], ldb);
		k++;
	    } else {

//2 x 2 diagonal block
//Interchange rows K+1 and -IPIV(K).
		kp = -ipiv[k];
		if (kp != k + 1) {
		    Rswap(nrhs, &B[k + 1 + ldb], ldb, &B[kp + ldb], ldb);
		}
//Multiply by inv(L(K)), where L(K) is the transformation
//stored in columns K and K+1 of A.

		if (k < n - 1) {
		    Rger(n - k - 1, nrhs, -One, &A[k + 2 + k * lda], 1, &B[k + ldb], ldb, &B[k + 2 + ldb], ldb);
		    Rger(n - k - 1, nrhs, -One, &A[k + 2 + (k + 1) * lda], 1, &B[k + 1 + ldb], ldb, &B[k + 2 + ldb], ldb);
		}
//Multiply by the inverse of the diagonal block.
		akm1k = A[k + 1 + k * lda];
		akm1 = A[k + k * lda] / akm1k;
		ak = A[k + 1 + (k + 1) * lda] / akm1k;
		denom = akm1 * ak - One;
		for (j = 0; j < nrhs; j++) {
		    bkm1 = B[k + j * ldb] / akm1k;
		    bk = B[k + 1 + j * ldb] / akm1k;
		    B[k + j * ldb] = (ak * bkm1 - bk) / denom;
		    B[k + 1 + j * ldb] = (akm1 * bk - bkm1) / denom;

		}
		k += 2;
	    }
	}

//Next solve L'*X = B, overwriting B with X.
//K is the main loop index, decreasing from N to 1 in steps of
//1 or 2, depending on the size of the diagonal blocks.
	k = n;
	while (1) {
//If K < 1, exit from loop.
	    if (k < 1)
		break;
	    if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Multiply by inv(L'(K)), where L(K) is the transformation
//stored in column K of A.
		if (k < n) {
		    Rgemv("Transpose", n - k, nrhs, -One, &B[k + 1 + ldb], ldb, &A[k + 1 + k * lda], 1, One, &B[k + ldb], ldb);
		}
//Interchange rows K and IPIV(K).
		kp = ipiv[k];
		if (kp != k) {
		    Rswap(nrhs, &B[k + ldb], ldb, &B[kp + ldb], ldb);
		}
		k--;
	    } else {

//2 x 2 diagonal block
//Multiply by inv(L'(K-1)), where L(K-1) is the transformation
//stored in columns K-1 and K of A.
		if (k < n) {
		    Rgemv("Transpose", n - k, nrhs, -One, &B[k + 1 + ldb], ldb, &A[k + 1 + k * lda], 1, One, &B[k + ldb], ldb);
		    Rgemv("Transpose", n - k, nrhs, -One, &B[k + 1 + ldb], ldb, &A[k + 1 + (k - 1) * lda], 1, One, &B[k - 1 + ldb], ldb);
		}
//Interchange rows K and -IPIV(K).
		kp = -ipiv[k];
		if (kp != k) {
		    Rswap(nrhs, &B[k + ldb], ldb, &B[kp + ldb], ldb);
		}
		k += -2;
	    }
	}
    }
    return;
}
