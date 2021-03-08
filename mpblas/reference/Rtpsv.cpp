/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rtpsv.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Rtpsv.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $

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
 Based on http://www.netlib.org/blas/dtpsv.f
Rtpsv solves one of the systems of equations
 A*x = b,   or   A'*x = b,
where b and x are n element vectors and A is an n by n unit, or
non-unit, upper or lower triangular matrix, supplied in packed form.
*/

#include <mpblas.h>

void Rtpsv(const char *uplo, const char *trans, const char *diag, INTEGER n, REAL * AP, REAL * x, INTEGER incx)
{
    INTEGER ix, j, jx, k, kx, kk, info, nounit;
    REAL Zero = 0.0;
    REAL temp;

//Test the input parameters.
    info = 0;
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
	Mxerbla("Rtpsv ", info);
	return;
    }
//quick return if possible.
    if (n == 0)
	return;

    nounit = Mlsame(diag, "N");

//set up the start point in x if the increment is not unity. this
//will be  ( n - 1 )*incx  too small for descending loops.
    kx = 0;
    if (incx <= 0)
	kx = -(n - 1) * incx;

//start the operations. in this version the elements of ap are
//accessed sequentially with one pass through ap.
    if (Mlsame(trans, "N")) {
//Form x := inv(A)*x.
	if (Mlsame(uplo, "U")) {
	    kk = (n * (n + 1)) / 2 - 1;
	    jx = kx + (n - 1) * incx;
	    for (j = n - 1; j >= 0; j--) {
		if (x[jx] != Zero) {
		    if (nounit)
			x[jx] = x[jx] / AP[kk];
		    temp = x[jx];
		    ix = jx;
		    for (k = kk - 1; k >= kk - j; k--) {
			ix = ix - incx;
			x[ix] = x[ix] - temp * AP[k];
		    }
		}
		jx = jx - incx;
		kk = kk - j - 1;
	    }
	} else {
	    kk = 0;
	    jx = kx;
	    for (j = 0; j < n; j++) {
		if (x[jx] != Zero) {
		    if (nounit)
			x[jx] = x[jx] / AP[kk];
		    temp = x[jx];
		    ix = jx;
		    for (k = kk + 1; k < kk + n - j; k++) {
			ix = ix + incx;
			x[ix] = x[ix] - temp * AP[k];
		    }
		}
		jx = jx + incx;
		kk = kk + n - j;
	    }
	}
    } else {
//Form x := inv(A')*x.
	if (Mlsame(uplo, "U")) {
	    kk = 0;
	    jx = kx;
	    for (j = 0; j < n; j++) {
		temp = x[jx];
		ix = kx;
		for (k = kk; k < kk + j; k++) {
		    temp = temp - AP[k] * x[ix];
		    ix = ix + incx;
		}
		if (nounit)
		    temp = temp / AP[kk + j];
		x[jx] = temp;
		jx = jx + incx;
		kk = kk + j + 1;
	    }
	} else {
	    kk = (n * (n + 1)) / 2 - 1;
	    kx = kx + (n - 1) * incx;
	    jx = kx;
	    for (j = n - 1; j >= 0; j--) {
		temp = x[jx];
		ix = kx;
		for (k = kk; k > kk - (n - j) + 1; k--) {
		    temp = temp - AP[k] * x[ix];
		    ix = ix - incx;
		}
		if (nounit)
		    temp = temp / AP[kk - (n - j) + 1];
		x[jx] = temp;
		jx = jx - incx;
		kk = kk - n + j;
	    }
	}
    }
    return;
}
