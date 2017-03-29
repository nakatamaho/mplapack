/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsptrf.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsptrf(const char *uplo, INTEGER n, REAL * ap, INTEGER * ipiv, INTEGER * info)
{
    INTEGER i, j, k;
    REAL t, r1, d11, d12, d21, d22;
    INTEGER kc, kk, kp;
    REAL wk;
    INTEGER kx, knc, kpc = 0, npp;
    REAL wkm1, wkp1;
    INTEGER imax = 0, jmax;
    REAL alpha;
    INTEGER kstep;
    INTEGER upper;
    REAL absakk;
    REAL colmax, rowmax;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	Mxerbla("Rsptrf", -(*info));
	return;
    }
//Initialize ALPHA for use in choosing pivot block size.
    alpha = (sqrt((REAL)17) + One) / (REAL)8;
    if (upper) {
//Factorize A as U*D*U' using the upper triangle of A
//K is the main loop index, decreasing from N to 1 in steps of
//1 or 2
	k = n;
	kc = (n - 1) * n / 2 + 1;
	while (1) {
	    knc = kc;

//If K < 1, exit from loop
	    if (k < 1)
		return;
	    kstep = 1;

//Determine rows and columns to be interchanged and whether
//a 1-by-1 or 2-by-2 pivot block will be used
	    absakk = abs(ap[kc + k - 1]);

//IMAX is the row-index of the largest off-diagonal element in
//column K, and COLMAX is its absolute value
	    if (k > 1) {
		imax = iRamax(k - 1, &ap[kc], 1);
		colmax = abs(ap[kc + imax - 1]);
	    } else {
		colmax = Zero;
	    }
	    if (max(absakk, colmax) == Zero) {
//Column K is zero: set INFO and continue
		if (*info == 0) {
		    *info = k;
		}
		kp = k;
	    } else {
		if (absakk >= alpha * colmax) {
//no interchange, use 1-by-1 pivot block
		    kp = k;
		} else {

//JMAX is the column-index of the largest off-diagonal
//element in row IMAX, and ROWMAX is its absolute value
		    rowmax = Zero;
		    jmax = imax;
		    kx = imax * (imax + 1) / 2 + imax;
		    for (j = imax + 1; j < k; j++) {
			if (abs(ap[kx]) > rowmax) {
			    rowmax = abs(ap[kx]);
			    jmax = j;
			}
			kx += j;
		    }
		    kpc = (imax - 1) * imax / 2 + 1;
		    if (imax > 1) {
			jmax = iRamax(imax - 1, &ap[kpc], 1);
			mtemp1 = rowmax, mtemp2 = abs(ap[kpc + jmax - 1]);
			rowmax = max(mtemp1, mtemp2);
		    }

		    if (absakk >= alpha * colmax * (colmax / rowmax)) {
//no interchange, use 1-by-1 pivot block
			kp = k;
		    } else if (abs(ap[kpc + imax - 1]) >= alpha * rowmax) {
//interchange rows and columns K and IMAX, use 1-by-1
//pivot block
			kp = imax;
		    } else {
//interchange rows and columns K-1 and IMAX, use 2-by-2
//pivot block
			kp = imax;
			kstep = 2;
		    }
		}
		kk = k - kstep + 1;
		if (kstep == 2) {
		    knc = knc - k + 1;
		}
		if (kp != kk) {
//Interchange rows and columns KK and KP in the leading
//submatrix A(1:k,1:k)

		    Rswap(kp - 1, &ap[knc], 1, &ap[kpc], 1);
		    kx = kpc + kp - 1;
		    for (j = kp + 1; j < kk - 1; j++) {
			kx = kx + j - 1;
			t = ap[knc + j - 1];
			ap[knc + j - 1] = ap[kx];
			ap[kx] = t;

		    }
		    t = ap[knc + kk - 1];
		    ap[knc + kk - 1] = ap[kpc + kp - 1];
		    ap[kpc + kp - 1] = t;
		    if (kstep == 2) {
			t = ap[kc + k - 2];
			ap[kc + k - 2] = ap[kc + kp - 1];
			ap[kc + kp - 1] = t;
		    }
		}
//Update the leading submatrix
		if (kstep == 1) {

//              1-by-1 pivot block D(k): column k now holds
//              W(k) = U(k)*D(k)
//              where U(k) is the k-th column of U
//              Perform a rank-1 update of A(1:k-1,1:k-1) as
//A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'

		    r1 = One / ap[kc + k - 1];
		    Rspr(uplo, k - 1, -r1, &ap[kc], 1, &ap[1]);

//Store U(k) in column k
		    Rscal(k - 1, r1, &ap[kc], 1);
		} else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )' */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )' */

		    if (k > 2) {

			d12 = ap[k - 1 + (k - 1) * k / 2];
			d22 = ap[k - 1 + (k - 2) * (k - 1) / 2] / d12;
			d11 = ap[k + (k - 1) * k / 2] / d12;
			t = One / (d11 * d22 - One);
			d12 = t / d12;

			for (j = k - 2; j >= 1; j--) {
			    wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2] - ap[j + (k - 1) * k / 2]);
			    wk = d12 * (d22 * ap[j + (k - 1) * k / 2] - ap[j + (k - 2) * (k - 1) / 2]);
			    for (i = j; i >= 1; i--) {
				ap[i + (j - 1) * j / 2] = ap[i + (j - 1) * j / 2] - ap[i + (k - 1) * k / 2] * wk - ap[i + (k - 2) * (k - 1) / 2] * wkm1;

			    }
			    ap[j + (k - 1) * k / 2] = wk;
			    ap[j + (k - 2) * (k - 1) / 2] = wkm1;

			}

		    }

		}
	    }

/*        Store details of the interchanges in IPIV */

	    if (kstep == 1) {
		ipiv[k] = kp;
	    } else {
		ipiv[k] = -kp;
		ipiv[k - 1] = -kp;
	    }

//Decrease K and return to the start of the main loop
	    k -= kstep;
	    kc = knc - k;
	}
    } else {

/*        Factorize A as L*D*L' using the lower triangle of A */
/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

	k = 0;
	kc = 1;
	npp = n * (n + 1) / 2;
	while (1) {
	    knc = kc;

//If K > N, exit from loop
	    if (k > n) {
		return;
	    }
	    kstep = 1;

//Determine rows and columns to be interchanged and whether
//a 1-by-1 or 2-by-2 pivot block will be used
	    absakk = abs(ap[kc]);

//IMAX is the row-index of the largest off-diagonal element in
//column K, and COLMAX is its absolute value
	    if (k < n) {
		imax = k + iRamax(n - k, &ap[kc + 1], 1);
		colmax = abs(ap[kc + imax - k]);
	    } else {
		colmax = Zero;
	    }

	    if (max(absakk, colmax) == Zero) {

//Column K is zero: set INFO and continue
		if (*info == 0) {
		    *info = k;
		}
		kp = k;
	    } else {
		if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

		    kp = k;
		} else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

		    rowmax = Zero;
		    kx = kc + imax - k;
		    for (j = k; j < imax - 1; j++) {
			if (abs(ap[kx]) > rowmax) {
			    rowmax = abs(ap[kx]);
			    jmax = j;
			}
			kx = kx + n - j;

		    }
		    kpc = npp - (n - imax + 1) * (n - imax + 2) / 2 + 1;
		    if (imax < n) {
			jmax = imax + iRamax(n - imax, &ap[kpc + 1], 1);
			mtemp1 = rowmax, mtemp2 = abs(ap[kpc + jmax - imax]);
			rowmax = max(mtemp1, mtemp2);
		    }

		    if (absakk >= alpha * colmax * (colmax / rowmax)) {
//no interchange, use 1-by-1 pivot block
			kp = k;
		    } else if (abs(ap[kpc]) >= alpha * rowmax) {

//interchange rows and columns K and IMAX, use 1-by-1
//pivot block
			kp = imax;
		    } else {

//interchange rows and columns K+1 and IMAX, use 2-by-2
//pivot block
			kp = imax;
			kstep = 2;
		    }
		}
		kk = k + kstep - 1;
		if (kstep == 2) {
		    knc = knc + n - k + 1;
		}
		if (kp != kk) {
//Interchange rows and columns KK and KP in the trailing
//submatrix A(k:n,k:n)
		    if (kp < n) {
			Rswap(n - kp, &ap[knc + kp - kk + 1], 1, &ap[kpc + 1], 1);
		    }
		    kx = knc + kp - kk;
		    for (j = kk + 1; j < kp - 1; j++) {
			kx = kx + n - j + 1;
			t = ap[knc + j - kk];
			ap[knc + j - kk] = ap[kx];
			ap[kx] = t;

		    }
		    t = ap[knc];
		    ap[knc] = ap[kpc];
		    ap[kpc] = t;
		    if (kstep == 2) {
			t = ap[kc + 1];
			ap[kc + 1] = ap[kc + kp - k];
			ap[kc + kp - k] = t;
		    }
		}
//Update the trailing submatrix
		if (kstep == 1) {

//1-by-1 pivot block D(k): column k now holds
//W(k) = L(k)*D(k)
//where L(k) is the k-th column of L
		    if (k < n) {

//Perform a rank-1 update of A(k+1:n,k+1:n) as
//A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
			r1 = One / ap[kc];
			Rspr(uplo, n - k, -r1, &ap[kc + 1], 1, &ap[kc + n - k + 1]);
//Store L(k) in column K
			Rscal(n - k, r1, &ap[kc + 1], 1);
		    }
		} else {

//2-by-2 pivot block D(k): columns K and K+1 now hold
//( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
//where L(k) and L(k+1) are the k-th and (k+1)-th columns
//of L

		    if (k < n - 1) {
//Perform a rank-2 update of A(k+2:n,k+2:n) as 
//A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
//= A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
			d21 = ap[k + 1 + (k - 1) * ((n << 1) - k) / 2];
			d11 = ap[k + 1 + k * ((n << 1) - k - 1) / 2] / d21;
			d22 = ap[k + (k - 1) * ((n << 1) - k) / 2] / d21;
			t = One / (d11 * d22 - One);
			d21 = t / d21;

			for (j = k + 2; j < n; j++) {
			    wk = d21 * (d11 * ap[j + (k - 1) * ((n << 1) - k) / 2] - ap[j + k * ((n << 1) - k - 1) / 2]);
			    wkp1 = d21 * (d22 * ap[j + k * ((n << 1) - k - 1) / 2] - ap[j + (k - 1) * ((n << 1) - k) / 2]);
			    for (i = j; i < n; i++) {
				ap[i + (j - 1) * ((n << 1) - j) / 2] = ap[i
									  + (j - 1) * ((n << 1) - j) / 2] - ap[i +
													       (k - 1) * ((n << 1) -
															  k) / 2] * wk - ap[i + k * ((n << 1) - k - 1) / 2] * wkp1;

			    }

			    ap[j + (k - 1) * ((n << 1) - k) / 2] = wk;
			    ap[j + k * ((n << 1) - k - 1) / 2] = wkp1;
			}
		    }
		}
	    }
//Store details of the interchanges in IPIV
	    if (kstep == 1) {
		ipiv[k] = kp;
	    } else {
		ipiv[k] = -kp;
		ipiv[k + 1] = -kp;
	    }

//Increase K and return to the start of the main loop
	    k += kstep;
	    kc = knc + n - k + 2;
	}

    }
    return;
}
