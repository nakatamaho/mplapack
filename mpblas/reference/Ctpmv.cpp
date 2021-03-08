/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Ctpmv.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Ctpmv.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $

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
http://www.netlib.org/blas/ctpmv.f
CTPMV  performs one of the matrix-vector operations
 x := A*x, or x := A'*x, or x := conjg( A' )*x,
where x is an n element vector and  A is an n by n unit, or non-unit,
upper or lower triangular matrix, supplied in packed form.
*/

#include <mpblas.h>

void Ctpmv(const char *uplo, const char *trans, const char *diag, INTEGER n, COMPLEX * AP, COMPLEX * x, INTEGER incx)
{
    INTEGER ix, j, jx, k, kx, kk, info = 0, noconj, nounit;
    REAL Zero = 0.0;
    COMPLEX temp;

//Test the input parameters.
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L"))
	info = 1;
    else if (!Mlsame(trans, "N") && !Mlsame(trans, "T") && !Mlsame(trans, "C"))
	info = 2;
    else if (!Mlsame(diag, "U") && !Mlsame(diag, "N"))
	info = 3;
    else if (n < 0)
	info = 4;
    else if (incx == 0)
	info = 7;
    if (info != 0) {
	Mxerbla("Ctpmv ", info);
	return;
    }
//quick return if possible.
    if (n == 0)
	return;
    noconj = Mlsame(trans, "T");
    nounit = Mlsame(diag, "N");
//set up the start point in x if the increment is not unity. this
//will be (n-1)*incx  too small for descending loops.
    kx = 0;
    if (incx <= 0)
	kx = (1 - n) * incx;

//start the operations. in this version the elements of ap are
//accessed sequentially with one pass through ap.
    if (Mlsame(trans, "N")) {
//Form x:= A*x.
	if (Mlsame(uplo, "U")) {
	    kk = 0;
	    jx = kx;
	    for (j = 0; j < n; j++) {
		if (x[jx] != Zero) {
		    temp = x[jx];
		    ix = kx;
		    for (k = kk; k < kk + j; k++) {
			x[ix] = x[ix] + temp * AP[k];
			ix = ix + incx;
		    }
		    if (nounit)
			x[jx] = x[jx] * AP[kk + j];
		}
		jx = jx + incx;
		kk = kk + j + 1;
	    }
	} else {
	    kk = (n * (n + 1)) / 2 - 1;
	    kx = kx + (n - 1) * incx;
	    jx = kx;
	    for (j = n - 1; j >= 0; j--) {
		if (x[jx] != Zero) {
		    temp = x[jx];
		    ix = kx;
		    for (k = kk; k > kk - (n - j) + 1; k--) {
			x[ix] = x[ix] + temp * AP[k];
			ix = ix - incx;
		    }
		    if (nounit)
			x[jx] = x[jx] * AP[kk - (n - j) + 1];
		}
		jx = jx - incx;
		kk = kk - (n - j);
	    }
	}
    } else {
	//Form x := A'*x or x := conjg(A')*x.
	if (Mlsame(uplo, "U")) {
	    kk = ((n * (n + 1)) / 2) - 1;
	    jx = kx + (n - 1) * incx;
	    for (j = n - 1; j >= 0; j--) {
		temp = x[jx];
		ix = jx;
		if (noconj) {
		    if (nounit)
			temp = temp * AP[kk];
		    for (k = kk - 1; k >= kk - j; k--) {
			ix = ix - incx;
			temp = temp + AP[k] * x[ix];
		    }
		} else {
		    if (nounit)
			temp = temp * conj(AP[kk]);
		    for (k = kk - 1; k >= kk - j; k--) {
			ix = ix - incx;
			temp = temp + conj(AP[k]) * x[ix];
		    }
		}
		x[jx] = temp;
		jx = jx - incx;
		kk = kk - j - 1;
	    }
	} else {
	    kk = 0;
	    jx = kx;
	    for (j = 0; j < n; j++) {
		temp = x[jx];
		ix = jx;
		if (noconj) {
		    if (nounit)
			temp = temp * AP[kk];
		    for (k = kk + 1; k < kk + n - j; k++) {
			ix = ix + incx;
			temp = temp + AP[k] * x[ix];
		    }
		} else {
		    if (nounit)
			temp = temp * conj(AP[kk]);
		    for (k = kk + 1; k < kk + n - j; k++) {
			ix = ix + incx;
			temp = temp + conj(AP[k]) * x[ix];
		    }
		}
		x[jx] = temp;
		jx = jx + incx;
		kk = kk + n - j;
	    }
	}
    }
    return;
}
