/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Ctrsv.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Ctrsv.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/ctrsv.f
Ctrsv solves one of the systems of equations
 A*x = b, or A'*x = b, or conjg( A' )*x = b,
where b and x are n element vectors and A is an n by n unit, or
non-unit, upper or lower triangular matrix.
*/

#include <mpblas.h>

void Ctrsv(const char *uplo, const char *trans, const char *diag, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * x, INTEGER incx)
{
    INTEGER i, ix, j, jx, kx, info = 0;
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
    else if (lda < max((INTEGER) 1, n))
	info = 6;
    else if (incx == 0)
	info = 8;
    if (info != 0) {
	Mxerbla("Ctrsv ", info);
	return;
    }
    //quick return if possible.
    if (n == 0)
	return;

    INTEGER noconj = Mlsame(trans, "T");
    INTEGER nounit = Mlsame(diag, "N");

//set up the start point in x if the increment is not unity. this
//will be (n-1)*incx too small for descending loops.
    if (incx <= 0)
	kx = (1 - n) * incx;
    else
	kx = 0;

//start the operations. in this version the elements of a are
//accessed sequentially with one pass through A.
    if (Mlsame(trans, "N")) {
//form  x := inv(A)*x.
	if (Mlsame(uplo, "U")) {
	    jx = kx + (n - 1) * incx;
	    for (j = n - 1; j >= 0; j--) {
		if (x[jx] != Zero) {
		    if (nounit)
			x[jx] = x[jx] / A[j + j * lda];
		    temp = x[jx];
		    ix = jx;
		    for (i = j - 1; i >= 0; i--) {
			ix = ix - incx;
			x[ix] = x[ix] - temp * A[i + j * lda];
		    }
		}
		jx = jx - incx;
	    }
	} else {
	    jx = kx;
	    for (j = 0; j < n; j++) {
		if (x[jx] != Zero) {
		    if (nounit)
			x[jx] = x[jx] / A[j + j * lda];
		    temp = x[jx];
		    ix = jx;
		    for (i = j + 1; i < n; i++) {
			ix = ix + incx;
			x[ix] = x[ix] - temp * A[i + j * lda];
		    }
		}
		jx = jx + incx;
	    }
	}
    } else {
//form x := inv(A')*x or x := inv(conjg(A'))*x.
	if (Mlsame(uplo, "U")) {
	    jx = kx;
	    for (j = 0; j < n; j++) {
		ix = kx;
		temp = x[jx];
		if (noconj) {
		    for (i = 0; i <= j - 1; i++) {
			temp = temp - A[i + j * lda] * x[ix];
			ix = ix + incx;
		    }
		    if (nounit)
			temp = temp / A[j + j * lda];
		} else {
		    for (i = 0; i <= j - 1; i++) {
			temp = temp - conj(A[i + j * lda]) * x[ix];
			ix = ix + incx;
		    }
		    if (nounit)
			temp = temp / conj(A[j + j * lda]);
		}
		x[jx] = temp;
		jx = jx + incx;
	    }
	} else {
	    kx = kx + (n - 1) * incx;
	    jx = kx;
	    for (j = n - 1; j >= 0; j--) {
		ix = kx;
		temp = x[jx];
		if (noconj) {
		    for (i = n - 1; i >= j + 1; i--) {
			temp = temp - A[i + j * lda] * x[ix];
			ix = ix - incx;
		    }
		    if (nounit)
			temp = temp / A[j + j * lda];
		} else {
		    for (i = n - 1; i >= j + 1; i--) {
			temp = temp - conj(A[i + j * lda]) * x[ix];
			ix = ix - incx;
		    }
		    if (nounit)
			temp = temp / conj(A[j + j * lda]);
		}
		x[jx] = temp;
		jx = jx - incx;
	    }
	}
    }
}
