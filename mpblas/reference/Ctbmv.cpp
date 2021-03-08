/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Ctbmv.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 *
 * $Id: Ctbmv.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $

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

/*
Based on /http://www.netlib.org/blas/ctbmv.f
Ctbmv performs one of the matrix-vector operations
x := A*x, or x := A'*x, or x := conjg(A')*x,
where x is an n element vector and  A is an n by n unit, or non-unit,
upper or lower triangular band matrix, with (k+1) diagonals.
*/

#include <mpblas.h>

void Ctbmv(const char *uplo, const char *trans, const char *diag, INTEGER n, INTEGER k, COMPLEX * A, INTEGER lda, COMPLEX * x,
	   INTEGER incx)
{
    INTEGER i, ix, j, jx, kx, l, info = 0;
    INTEGER noconj;
    INTEGER nounit;
    COMPLEX temp;
    REAL Zero = 0.0;

//Test the input parameters.
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L"))
	info = 1;
    else if (!Mlsame(trans, "N") && !Mlsame(trans, "T")
	     && !Mlsame(trans, "C"))
	info = 2;
    else if (!Mlsame(diag, "U") && !Mlsame(diag, "N"))
	info = 3;
    else if (n < 0)
	info = 4;
    else if (k < 0)
	info = 5;
    else if (lda < (k + 1))
	info = 7;
    else if (incx == 0)
	info = 9;
    if (info != 0) {
	Mxerbla("Ctbmv ", info);
	return;
    }
//quick return if possible.
    if (n == 0)
	return;
    noconj = Mlsame(trans, "T");
    nounit = Mlsame(diag, "N");

//set up the start point in x if the increment is not unity. this
//will be (n-1)*incx too small for descending loops.
    kx = 0;
    if (incx <= 0)
	kx = -(n - 1) * incx;

//start the operations. in this version the elements of A are
//accessed sequentially with one pass through A.
    if (Mlsame(trans, "N")) {
//Form x := A*x.
	if (Mlsame(uplo, "U")) {
	    jx = kx;
	    for (j = 0; j < n; j++) {
		if (x[jx] != Zero) {
		    temp = x[jx];
		    ix = kx;
		    for (i = max((INTEGER) 0, j - k); i < j; i++) {
			x[ix] = x[ix] + temp * A[(k - j + i) + j * lda];
			ix = ix + incx;
		    }
		    if (nounit)
			x[jx] = x[jx] * A[k + j * lda];
		}
		jx = jx + incx;
		if (j >= k)
		    kx = kx + incx;
	    }
	} else {
	    kx = kx + (n - 1) * incx;
	    jx = kx;
	    for (j = n - 1; j >= 0; j--) {
		if (x[jx] != Zero) {
		    temp = x[jx];
		    ix = kx;
		    l = 1 - j;
		    for (i = min(n, j + k + 1) - 1; i >= j + 1; i--) {
			x[ix] = x[ix] + temp * A[(i - j) + j * lda];
			ix = ix - incx;
		    }
		    if (nounit)
			x[jx] = x[jx] * A[j * lda];
		}
		jx = jx - incx;
		if ((n - j) > k)
		    kx = kx - incx;
	    }
	}
    } else {
//form x := A'*x or x := conjg(A')*x.
	if (Mlsame(uplo, "U")) {
	    kx = kx + (n - 1) * incx;
	    jx = kx;
	    for (j = n - 1; j >= 0; j--) {
		temp = x[jx];
		kx = kx - incx;
		ix = kx;
		if (noconj) {
		    if (nounit)
			temp = temp * A[k + j * lda];
		    for (i = j - 1; i >= max((INTEGER) 0, j - k); i--) {
			temp = temp + A[(k - j + i) + j * lda] * x[ix];
			ix = ix - incx;
		    }
		} else {
		    if (nounit)
			temp = temp * conj(A[k + j * lda]);
		    for (i = j - 1; i >= max((INTEGER) 0, j - k); i--) {
			temp = temp + conj(A[(k - j + i) + j * lda]) * x[ix];
			ix = ix - incx;
		    }
		}
		x[jx] = temp;
		jx = jx - incx;
	    }
	} else {
	    jx = kx;
	    for (j = 0; j < n; j++) {
		temp = x[jx];
		kx = kx + incx;
		ix = kx;
		if (noconj) {
		    if (nounit)
			temp = temp * A[j * lda];
		    for (i = j + 1; i < min(n, j + k + 1); i++) {
			temp = temp + A[(i - j) + j * lda] * x[ix];
			ix = ix + incx;
		    }
		} else {
		    if (nounit)
			temp = temp * conj(A[j * lda]);
		    for (i = j + 1; i < min(n, j + k + 1); i++) {
			temp = temp + conj(A[(i - j) + j * lda]) * x[ix];
			ix = ix + incx;
		    }
		}
		x[jx] = temp;
		jx = jx + incx;
	    }
	}
    }
    return;
}
