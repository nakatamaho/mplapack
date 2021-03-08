/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Csytri.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Csytri(const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * ipiv, COMPLEX * work, INTEGER * info)
{
    COMPLEX d;
    INTEGER k;
    COMPLEX t, ak;
    INTEGER kp;
    COMPLEX akp1, temp, akkp1;
    INTEGER kstep;
    INTEGER upper;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Csytri", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Check that the diagonal matrix D is nonsingular.
    if (upper) {
//Upper triangular storage: examine D from bottom to top
	for (*info = n; *info >= 1; --(*info)) {
	    if (ipiv[*info] > 0 && A[*info + *info * lda] == Zero) {
		return;
	    }
	}
    } else {
//Lower triangular storage: examine D from top to bottom.
	for (*info = 1; *info <= n; ++(*info)) {
	    if (ipiv[*info] > 0 && A[*info + *info * lda] == Zero) {
		return;
	    }
	}
    }
    *info = 0;
    if (upper) {
//Compute inv(A) from the factorization A = U*D*U'.
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2, depending on the size of the diagonal blocks.
	k = 0;
      L30:
//If K > N, exit from loop.
	if (k > n) {
	    goto L40;
	}
	if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Invert the diagonal block.
	    A[k + k * lda] = One / A[k + k * lda];
//Compute column K of the inverse.
	    if (k > 1) {
		Ccopy(k - 1, &A[k * lda], 1, &work[0], 1);
		Csymv(uplo, k - 1, (COMPLEX) - One, &A[0], lda, &work[0], 1, (COMPLEX) Zero, &A[k * lda], 1);
		A[k + k * lda] = A[k + k * lda] - Cdotu(k - 1, &work[0], 1, &A[k * lda], 1);
	    }
	    kstep = 1;
	} else {
//2 x 2 diagonal block
//Invert the diagonal block
	    t = A[k + (k + 1) * lda];
	    ak = A[k + k * lda] / t;
	    akp1 = A[k + 1 + (k + 1) * lda] / t;
	    akkp1 = A[k + (k + 1) * lda] / t;
	    d = t * (ak * akp1 - One);
	    A[k + k * lda] = akp1 / d;
	    A[k + 1 + (k + 1) * lda] = ak / d;
	    A[k + (k + 1) * lda] = A[k + (k + 1) * lda] - akkp1 / d;
//Compute columns K and K+1 of the inverse.
	    if (k > 1) {
		Ccopy(k - 1, &A[k * lda], 1, &work[0], 1);
		Csymv(uplo, k - 1, (COMPLEX) - One, &A[0], lda, &work[0], 1, (COMPLEX) Zero, &A[k * lda], 1);
		A[k + k * lda] = A[k + k * lda] - Cdotu(k - 1, &work[0], 1, &A[k * lda], 1);
		A[k + (k + 1) * lda] = A[k + (k + 1) * lda] - Cdotu(k - 1, &A[k * lda], 1, &A[(k + 1) * lda], 1);
		Ccopy(k - 1, &A[(k + 1) * lda], 1, &work[0], 1);
		Csymv(uplo, k - 1, (COMPLEX) - One, &A[0], lda, &work[0], 1, (COMPLEX) Zero, &A[(k + 1) * lda], 1);
		A[k + 1 + (k + 1) * lda] = A[k + 1 + (k + 1) * lda] - Cdotu(k - 1, &work[0], 1, &A[(k + 1) * lda], 1);
	    }
	    kstep = 2;
	}
	kp = abs(ipiv[k]);
	if (kp != k) {
//Interchange rows and columns K and KP in the leading
//submatrix A(1:k+1,1:k+1)
	    Cswap(kp - 1, &A[k * lda], 1, &A[kp * lda], 1);
	    Cswap(k - kp - 1, &A[kp + 1 + k * lda], 1, &A[kp + (kp + 1) * lda], lda);
	    temp = A[k + k * lda];
	    A[k + k * lda] = A[kp + kp * lda];
	    A[kp + kp * lda] = temp;
	    if (kstep == 2) {
		temp = A[k + (k + 1) * lda];
		A[k + (k + 1) * lda] = A[kp + (k + 1) * lda];
		A[kp + (k + 1) * lda] = temp;
	    }
	}
	k = k + kstep;
	goto L30;
      L40:
	;
    } else {
//Compute inv(A) from the factorization A = L*D*L'.
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2, depending on the size of the diagonal blocks.
	k = n;
      L50:
//If K < 1, exit from loop.
	if (k < 1) {
	    goto L60;
	}
	if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Invert the diagonal block.
	    A[k + k * lda] = One / A[k + k * lda];
//Compute column K of the inverse.
	    if (k < n) {
		Ccopy(n - k, &A[k + 1 + k * lda], 1, &work[0], 1);
		Csymv(uplo, n - k, (COMPLEX) - One, &A[k + 1 + (k + 1) * lda], lda, &work[0], 1, (COMPLEX) Zero, &A[k + 1 + k * lda], 1);
		A[k + k * lda] = A[k + k * lda] - Cdotu(n - k, &work[0], 1, &A[k + 1 + k * lda], 1);
	    }
	    kstep = 1;
	} else {
//2 x 2 diagonal block
//Invert the diagonal block.
	    t = A[k + (k - 1) * lda];
	    ak = A[k - 1 + (k - 1) * lda] / t;
	    akp1 = A[k + k * lda] / t;
	    akkp1 = A[k + (k - 1) * lda] / t;
	    d = t * (ak * akp1 - One);
	    A[k - 1 + (k - 1) * lda] = akp1 / d;
	    A[k + k * lda] = ak / d;
	    A[k + (k - 1) * lda] = -akkp1 / d;
//Compute columns K-1 and K of the inverse.
	    if (k < n) {
		Ccopy(n - k, &A[k + 1 + k * lda], 1, &work[0], 1);
		Csymv(uplo, n - k, (COMPLEX) - One, &A[k + 1 + (k + 1) * lda], lda, &work[0], 1, (COMPLEX) Zero, &A[k + 1 + k * lda], 1);
		A[k + k * lda] = A[k + k * lda] - Cdotu(n - k, &work[0], 1, &A[k + 1 + k * lda], 1);
		A[k + (k - 1) * lda] = A[k + (k - 1) * lda] - Cdotu(n - k, &A[k + 1 + k * lda], 1, &A[k + 1 + (k - 1) * lda], 1);
		Ccopy(n - k, &A[k + 1 + (k - 1) * lda], 1, &work[0], 1);
		Csymv(uplo, n - k, (COMPLEX) - One, &A[k + 1 + (k + 1) * lda], lda, &work[0], 1, (COMPLEX) Zero, &A[k + 1 + (k - 1) * lda], 1);
		A[k - 1 + (k - 1) * lda] = A[k - 1 + (k - 1) * lda] - Cdotu(n - k, &work[0], 1, &A[k + 1 + (k - 1) * lda], 1);
	    }
	    kstep = 2;
	}
	kp = abs(ipiv[k]);
	if (kp != k) {
//Interchange rows and columns K and KP in the trailing
//submatrix A(k-1:n,k-1:n)
	    if (kp < n) {
		Cswap(n - kp, &A[kp + 1 + k * lda], 1, &A[kp + 1 + kp * lda], 1);
	    }
	    Cswap(kp - k - 1, &A[k + 1 + k * lda], 1, &A[kp + (k + 1) * lda], lda);
	    temp = A[k + k * lda];
	    A[k + k * lda] = A[kp + kp * lda];
	    A[kp + kp * lda] = temp;
	    if (kstep == 2) {
		temp = A[k + (k - 1) * lda];
		A[k + (k - 1) * lda] = A[kp + (k - 1) * lda];
		A[kp + (k - 1) * lda] = temp;
	    }
	}
	k = k - kstep;
	goto L50;
      L60:
	;
    }
    return;
}
