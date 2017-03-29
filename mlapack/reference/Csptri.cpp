/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Csptri.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Csptri(const char *uplo, INTEGER n, COMPLEX * ap, INTEGER * ipiv, COMPLEX * work, INTEGER * info)
{
    COMPLEX d;
    INTEGER j, k;
    COMPLEX t, ak;
    INTEGER kc, kp, kx, kpc, npp;
    COMPLEX akp1, temp, akkp1;
    INTEGER kstep;
    INTEGER upper;
    INTEGER kcnext;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	Mxerbla("Csptri", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Check that the diagonal matrix D is nonsingular.
    if (upper) {
//Upper triangular storage: examine D from bottom to top
	kp = n * (n + 1) / 2;
	for (*info = n; *info >= 1; --(*info)) {
	    if (ipiv[*info] > 0 && ap[kp] == Zero) {
		return;
	    }
	    kp = kp - *info;
	}
    } else {
//Lower triangular storage: examine D from top to bottom.
	kp = 1;
	for (*info = 1; *info <= n; ++(*info)) {
	    if (ipiv[*info] > 0 && ap[kp] == Zero) {
		return;
	    }
	    kp = kp + n - *info + 1;
	}
    }
    *info = 0;
    if (upper) {
//Compute inv(A) from the factorization A = U*D*U'.
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2, depending on the size of the diagonal blocks.
	k = 0;
	kc = 1;
      L30:
//If K > N, exit from loop.
	if (k > n) {
	    goto L50;
	}
	kcnext = kc + k;
	if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Invert the diagonal block.
	    ap[kc + k - 1] = One / ap[kc + k - 1];
//Compute column K of the inverse.
	    if (k > 1) {
		Ccopy(k - 1, &ap[kc], 1, &work[0], 1);
		Cspmv(uplo, k - 1, (COMPLEX) - One, &ap[1], &work[0], 1, (COMPLEX) Zero, &ap[kc], 1);
		ap[kc + k - 1] = ap[kc + k - 1] - Cdotu(k - 1, &work[0], 1, &ap[kc], 1);
	    }
	    kstep = 1;
	} else {
//2 x 2 diagonal block
//Invert the diagonal block.
	    t = ap[kcnext + k - 1];
	    ak = ap[kc + k - 1] / t;
	    akp1 = ap[kcnext + k] / t;
	    akkp1 = ap[kcnext + k - 1] / t;
	    d = t * (ak * akp1 - One);
	    ap[kc + k - 1] = akp1 / d;
	    ap[kcnext + k] = ak / d;
	    ap[kcnext + k - 1] = -akkp1 / d;
//Compute columns K and K+1 of the inverse.
	    if (k > 1) {
		Ccopy(k - 1, &ap[kc], 1, &work[0], 1);
		Cspmv(uplo, k - 1, (COMPLEX) - One, &ap[1], &work[0], 1, Zero, &ap[kc], 1);
		ap[kc + k - 1] = ap[kc + k - 1] - Cdotu(k - 1, &work[0], 1, &ap[kc], 1);
		ap[kcnext + k - 1] = ap[kcnext + k - 1] - Cdotu(k - 1, &ap[kc], 1, &ap[kcnext], 1);
		Ccopy(k - 1, &ap[kcnext], 1, &work[0], 1);
		Cspmv(uplo, k - 1, (COMPLEX) - One, &ap[1], &work[0], 1, Zero, &ap[kcnext], 1);
		ap[kcnext + k] = ap[kcnext + k] - Cdotu(k - 1, &work[0], 1, &ap[kcnext], 1);
	    }
	    kstep = 2;
	    kcnext = kcnext + k + 1;
	}
	kp = abs(ipiv[k]);
	if (kp != k) {
//Interchange rows and columns K and KP in the leading
//submatrix A(1:k+1,1:k+1)
	    kpc = (kp - 1) * kp / 2 + 1;
	    Cswap(kp - 1, &ap[kc], 1, &ap[kpc], 1);
	    kx = kpc + kp - 1;
	    for (j = kp + 1; j <= k - 1; j++) {
		kx = kx + j - 1;
		temp = ap[kc + j - 1];
		ap[kc + j - 1] = ap[kx];
		ap[kx] = temp;
	    }
	    temp = ap[kc + k - 1];
	    ap[kc + k - 1] = ap[kpc + kp - 1];
	    ap[kpc + kp - 1] = temp;
	    if (kstep == 2) {
		temp = ap[kc + k + k - 1];
		ap[kc + k + k - 1] = ap[kc + k + kp - 1];
		ap[kc + k + kp - 1] = temp;
	    }
	}
	k = k + kstep;
	kc = kcnext;
	goto L30;
      L50:
	;
    } else {
//Compute inv(A) from the factorization A = L*D*L'.
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2, depending on the size of the diagonal blocks.
	npp = n * (n + 1) / 2;
	k = n;
	kc = npp;
      L60:
//If K < 1, exit from loop.
	if (k < 1) {
	    goto L80;
	}
	kcnext = kc - (n - k + 2);
	if (ipiv[k] > 0) {
//1 x 1 diagonal block
//Invert the diagonal block.
	    ap[kc] = One / ap[kc];
//Compute column K of the inverse.
	    if (k < n) {
		Ccopy(n - k, &ap[kc + 1], 1, &work[0], 1);
		Cspmv(uplo, n - k, (COMPLEX) - One, &ap[kc + n - k + 1], &work[0], 1, (COMPLEX) Zero, &ap[kc + 1], 1);
		ap[kc] = ap[kc] - Cdotu(n - k, &work[0], 1, &ap[kc + 1], 1);
	    }
	    kstep = 1;
	} else {
//2 x 2 diagonal block
//Invert the diagonal block.
	    t = ap[kcnext + 1];
	    ak = ap[kcnext] / t;
	    akp1 = ap[kc] / t;
	    akkp1 = ap[kcnext + 1] / t;
	    d = t * (ak * akp1 - One);
	    ap[kcnext] = akp1 / d;
	    ap[kc] = ak / d;
	    ap[kcnext + 1] = -akkp1 / d;
//Compute columns K-1 and K of the inverse.
	    if (k < n) {
		Ccopy(n - k, &ap[kc + 1], 1, &work[0], 1);
		Cspmv(uplo, n - k, (COMPLEX) - One, &ap[kc + (n - k + 1)], &work[0], 1, Zero, &ap[kc + 1], 1);
		ap[kc] = ap[kc] - Cdotu(n - k, &work[0], 1, &ap[kc + 1], 1);
		ap[kcnext + 1] = ap[kcnext + 1] - Cdotu(n - k, &ap[kc + 1], 1, &ap[kcnext + 2], 1);
		Ccopy(n - k, &ap[kcnext + 2], 1, &work[0], 1);
		Cspmv(uplo, n - k, (COMPLEX) - One, &ap[kc + (n - k + 1)], &work[0], 1, Zero, &ap[kcnext + 2], 1);
		ap[kcnext] = ap[kcnext] - Cdotu(n - k, &work[0], 1, &ap[kcnext + 2], 1);
	    }
	    kstep = 2;
	    kcnext = kcnext - n - k + 3;
	}
	kp = abs(ipiv[k]);
	if (kp != k) {
//Interchange rows and columns K and KP in the trailing
//submatrix A(k-1:n,k-1:n)
	    kpc = npp - (n - kp + 1) * (n - kp + 2) / 2 + 1;
	    if (kp < n) {
		Cswap(n - kp, &ap[kc + kp - k + 1], 1, &ap[kpc + 1], 1);
	    }
	    kx = kc + kp - k;
	    for (j = k + 1; j <= kp - 1; j++) {
		kx = kx + n - j + 1;
		temp = ap[kc + j - k];
		ap[kc + j - k] = ap[kx];
		ap[kx] = temp;

	    }
	    temp = ap[kc];
	    ap[kc] = ap[kpc];
	    ap[kpc] = temp;
	    if (kstep == 2) {
		temp = ap[kc - n + k - 1];
		ap[kc - n + k - 1] = ap[kc - n + kp - 1];
		ap[kc - n + kp - 1] = temp;
	    }
	}
	k = k - kstep;
	kc = kcnext;
	goto L60;
      L80:
	;
    }
    return;
}
