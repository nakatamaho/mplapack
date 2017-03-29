/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chetf2.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chetf2(const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * ipiv, INTEGER * info)
{
    REAL d;
    INTEGER i, j, k;
    COMPLEX t;
    REAL r1, d11;
    COMPLEX d12;
    REAL d22;
    COMPLEX d21;
    INTEGER kk, kp;
    COMPLEX wk;
    REAL tt;
    COMPLEX wkm1, wkp1;
    INTEGER imax = 0, jmax;
    REAL alpha;
    INTEGER kstep;
    INTEGER upper;
    REAL absakk;
    REAL colmax;
    REAL rowmax;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2 = 0.0;

//Test the input parameters.
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
	Mxerbla("Chetf2", -(*info));
	return;
    }
//Initialize ALPHA for use in choosing pivot block size.
    alpha = (sqrt((REAL)17) + One) / (REAL)8;
    if (upper) {
//Factorize A as U*D*U' using the upper triangle of A
//K is the main loop index, decreasing from N to 1 in steps of
//1 or 2
	k = n;
      L10:
//If K < 1, exit from loop
	if (k < 1) {
	    goto L90;
	}
	kstep = 1;
//Determine rows and columns to be interchanged and whether
//a 1-by-1 or 2-by-2 pivot block will be used
	absakk = abs(A[k + k * lda].real());
//IMAX is the row-index of the largest off-diagonal element in
//column K, and COLMAX is its absolute value
	if (k > 1) {
	    imax = iCamax(k - 1, &A[k * lda], 1);
	    colmax = abs(A[imax + k * lda].real()) + abs(A[imax + k * lda].imag());
	} else {
	    colmax = Zero;
	}

	if (max(absakk, colmax) == Zero) {
//Column K is zero or contains a NaN: set INFO and continue
	    if (*info == 0) {
		*info = k;
	    }
	    kp = k;
	    A[k + k * lda] = A[k + k * lda].real();
	} else {
	    if (absakk >= alpha * colmax) {
//no interchange, use 1-by-1 pivot block
		kp = k;
	    } else {
//JMAX is the column-index of the largest off-diagonal
//element in row IMAX, and ROWMAX is its absolute value
		jmax = imax + iCamax(k - imax, &A[imax + (imax + 1) * lda], lda);
		rowmax = abs(A[imax + jmax * lda].real()) + abs(A[imax + jmax * lda].imag());
		if (imax > 1) {
		    jmax = iCamax(imax - 1, &A[imax * lda], 1);
		    mtemp1 = rowmax, mtemp2 = abs(A[jmax + imax * lda].real()) + abs(A[jmax + imax * lda].imag());
		    rowmax = max(mtemp1, mtemp2);
		}
		if (absakk >= alpha * colmax * (colmax / rowmax)) {
//no interchange, use 1-by-1 pivot block
		    kp = k;
		} else if (abs(A[imax + imax * lda].real()) >= alpha * rowmax) {
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
	    if (kp != kk) {
//Interchange rows and columns KK and KP in the leading
//submatrix A(1:k,1:k)
		Cswap(kp - 1, &A[kk * lda], 1, &A[kp * lda], 1);
		for (j = kp + 1; j <= kk - 1; j++) {
		    t = conj(A[j + kk * lda]);
		    A[j + kk * lda] = conj(A[kp + j * lda]);
		    A[kp + j * lda] = t;
		}
		A[kp + kk * lda] = conj(A[kp + kk * lda]);
		r1 = A[kk + kk * lda].real();
		A[kk + kk * lda] = A[kp + kp * lda].real();
		A[kp + kp * lda] = r1;
		if (kstep == 2) {
		    A[k + k * lda] = A[k + k * lda].real();
		    t = A[k - 1 + k * lda];
		    A[k - 1 + k * lda] = A[kp + k * lda];
		    A[kp + k * lda] = t;
		}
	    } else {
		A[k + k * lda] = A[k + k * lda].real();
		if (kstep == 2) {
		    A[k - 1 + (k - 1) * lda] = A[k - 1 + (k - 1) * lda].real();
		}
	    }
//Update the leading submatrix
	    if (kstep == 1) {
//1-by-1 pivot block D(k): column k now holds
//W(k) = U(k)*D(k)
//where U(k) is the k-th column of U
//Perform a rank-1 update of A(1:k-1,1:k-1) as
//A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
		Cher(uplo, k - 1, -One / A[k + k * lda].real(), &A[k * lda], 1, &A[0], lda);
//Store U(k) in column k
		CRscal(k - 1, -One, &A[k * lda], 1);
	    } else {
//2-by-2 pivot block D(k): columns k and k-1 now hold
//( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
//where U(k) and U(k-1) are the k-th and (k-1)-th columns
//of U
//Perform a rank-2 update of A(1:k-2,1:k-2) as
//A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
//   = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
		if (k > 2) {
		    d = Rlapy2(A[k - 1 + k * lda].real(), A[k - 1 + k * lda].imag());
		    d22 = A[k - 1 + (k - 1) * lda].real() / d;
		    d11 = A[k + k * lda].real() / d;
		    tt = One / (d11 * d22 - One);
		    d12 = A[k - 1 + k * lda] / d;
		    d = tt / d;

		    for (j = k - 2; j >= 1; j--) {
			wkm1 = d * (d11 * A[j + (k - 1) * lda] - conj(d12) * A[j + k * lda]);
			wk = d * (d22 * A[j + k * lda] - d12 * A[j + (k - 1) * lda]);
			for (i = j; i >= 1; i--) {
			    A[i + j * lda] = A[i + j * lda] - A[i + k * lda] * conj(wk) - A[i + (k - 1) * lda] * conj(wkm1);
			}
			A[j + k * lda] = wk;
			A[j + (k - 1) * lda] = wkm1;
			A[j + j * lda] = A[j + j * lda].real();
		    }
		}
	    }
	}
//Store details of the interchanges in IPIV
	if (kstep == 1) {
	    ipiv[k] = kp;
	} else {
	    ipiv[k] = -kp;
	    ipiv[k - 1] = -kp;
	}
//Decrease K and return to the start of the main loop
	k = k - kstep;
	goto L10;
    } else {
//Factorize A as L*D*L' using the lower triangle of A
//K is the main loop index, increasing from 1 to N in steps of
//1 or 2
	k = 0;
      L50:
//If K > N, exit from loop
	if (k > n) {
	    goto L90;
	}
	kstep = 1;
//Determine rows and columns to be interchanged and whether
//a 1-by-1 or 2-by-2 pivot block will be used
	absakk = abs(A[k + k * lda].real());

//IMAX is the row-index of the largest off-diagonal element in
//column K, and COLMAX is its absolute value
	if (k < n) {
	    imax = k + iCamax(n - k, &A[k + 1 + k * lda], 1);
	    colmax = abs(A[imax + k * lda].real()) + abs(A[imax + k * lda].imag());
	} else {
	    colmax = Zero;
	}
	if (max(absakk, colmax) == Zero) {
//Column K is zero or contains a NaN: set INFO and continue
	    if (*info == 0) {
		*info = k;
	    }
	    kp = k;
	    A[k + k * lda] = A[k + k * lda].real();
	} else {
	    if (absakk >= alpha * colmax) {
//no interchange, use 1-by-1 pivot block
		kp = k;
	    } else {
//JMAX is the column-index of the largest off-diagonal
//element in row IMAX, and ROWMAX is its absolute value
		jmax = k - 1 + iCamax(imax - k, &A[imax + k * lda], lda);
		rowmax = abs(A[imax + jmax * lda].real()) + abs(A[imax + jmax * lda].imag());
		if (imax < n) {
		    jmax = imax + iCamax(n - imax, &A[imax + 1 + imax * lda], 1);
		    mtemp1 = rowmax, mtemp1 = abs(A[jmax + imax * lda].real()) + abs(A[jmax + imax * lda].imag());
		    rowmax = max(mtemp1, mtemp2);
		}
		if (absakk >= alpha * colmax * (colmax / rowmax)) {
//no interchange, use 1-by-1 pivot block
		    kp = k;
		} else if (abs(A[imax + imax * lda].real()) >= alpha * rowmax) {
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
	    if (kp != kk) {
//Interchange rows and columns KK and KP in the trailing
//submatrix A(k:n,k:n)
		if (kp < n) {
		    Cswap(n - kp, &A[kp + 1 + kk * lda], 1, &A[kp + 1 + kp * lda], 1);
		}
		for (j = kk + 1; j <= kp - 1; j++) {
		    t = conj(A[j + kk * lda]);
		    A[j + kk * lda] = conj(A[kp + j * lda]);
		    A[kp + j * lda] = t;
		}
		A[kp + kk * lda] = conj(A[kp + kk * lda]);
		r1 = A[kk + kk * lda].real();
		A[kk + kk * lda] = A[kp + kp * lda].real();
		A[kp + kp * lda] = r1;
		if (kstep == 2) {
		    A[k + k * lda] = A[k + k * lda].real();
		    t = A[k + 1 + k * lda];
		    A[k + 1 + k * lda] = A[kp + k * lda];
		    A[kp + k * lda] = t;
		}
	    } else {
		A[k + k * lda] = A[k + k * lda].real();
		if (kstep == 2) {
		    A[k + 1 + (k + 1) * lda] = A[k + 1 + (k + 1) * lda].real();
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
		    r1 = One / (REAL) A[k + k * lda].real();
		    Cher(uplo, n - k, -r1, &A[k + 1 + k * lda], 1, &A[k + 1 + (k + 1) * lda], lda);
//Store L(k) in column K
		    CRscal(n - k, r1, &A[k + 1 + k * lda], 1);
		}
	    } else {
//2-by-2 pivot block D(k)
		if (k < n - 1) {
//Perform a rank-2 update of A(k+2:n,k+2:n) as
//A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
//   = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
//  where L(k) and L(k+1) are the k-th and (k+1)-th
//  columns of L
		    d = Rlapy2(A[k + 1 + k * lda].real(), A[k + 1 + k * lda].imag());
		    d11 = A[k + 1 + (k + 1) * lda].real() / d;
		    d22 = A[k + k * lda].real() / d;
		    tt = One / (d11 * d22 - One);
		    d21 = A[k + 1 + k * lda] / d;
		    d = tt / d;
		    for (j = k + 2; j <= n; j++) {
			wk = d * (d11 * A[j + k * lda] - d21 * A[j + (k + 1) * lda]);
			wkp1 = d * (d22 * A[j + (k + 1) * lda] - conj(d21) * A[j + k * lda]);
			for (i = j; i <= n; i++) {
			    A[i + j * lda] = A[i + j * lda] - A[i + k * lda] * conj(wk) - A[i + (k + 1) * lda] * conj(wkp1);

			}
			A[j + k * lda] = wk;
			A[j + (k + 1) * lda] = wkp1;
			A[j + j * lda] = A[j + j * lda].real();

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
	goto L50;

    }
  L90:
    return;
}
