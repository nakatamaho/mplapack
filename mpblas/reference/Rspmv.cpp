/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rspmv.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Rspmv.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/dspmv.f
Rspmv performs the matrix-vector operation
y := alpha*A*x + beta*y,
where alpha and beta are scalars, x and y are n element vectors and
A is an n by n symmetric matrix, supplied in packed form.
*/

#include <mpblas.h>

void Rspmv(const char *uplo, INTEGER n, REAL alpha, REAL * AP, REAL * x, INTEGER incx, REAL beta, REAL * y, INTEGER incy)
{
    INTEGER i, ix, iy, j, jx, jy, k, kx, ky, kk, info = 0;
    REAL Zero = 0.0, One = 1.0;
    REAL temp1, temp2;

//Test the input parameters.
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L"))
	info = 1;
    else if (n < 0)
	info = 2;
    else if (incx == 0)
	info = 6;
    else if (incy == 0)
	info = 9;

    if (info != 0) {
	Mxerbla("Rspmv ", info);
	return;
    }
//quick return if possible.
    if ((n == 0) || ((alpha == Zero) && (beta == One)))
	return;

//set up the start points in  x  and  y.
    if (incx > 0)
	kx = 0;
    else
	kx = -(n - 1) * incx;
    if (incy > 0)
	ky = 0;
    else
	ky = -(n - 1) * incy;

//Start the operations. In this version the elements of the array AP
//are accessed sequentially with one pass through AP.
//First form  y := beta*y.
    if (beta != One) {
	iy = ky;
	if (beta == Zero) {
	    for (i = 0; i < n; i++) {
		y[iy] = Zero;
		iy = iy + incy;
	    }
	} else {
	    for (i = 0; i < n; i++) {
		y[iy] = beta * y[iy];
		iy = iy + incy;
	    }
	}
    }
    if (alpha == Zero)
	return;

    kk = 0;
    if (Mlsame(uplo, "U")) {
//Form y when AP contains the upper triangle.
	jx = kx;
	jy = ky;
	for (j = 0; j < n; j++) {
	    temp1 = alpha * x[jx];
	    temp2 = Zero;
	    ix = kx;
	    iy = ky;
	    for (k = kk; k < kk + j; k++) {
		y[iy] = y[iy] + temp1 * AP[k];
		temp2 = temp2 + AP[k] * x[ix];
		ix = ix + incx;
		iy = iy + incy;
	    }
	    y[jy] = y[jy] + temp1 * AP[kk + j] + alpha * temp2;
	    jx = jx + incx;
	    jy = jy + incy;
	    kk = kk + j + 1;
	}
    } else {
//Form y when AP contains the lower triangle.
	jx = kx;
	jy = ky;
	for (j = 0; j < n; j++) {
	    temp1 = alpha * x[jx];
	    temp2 = Zero;
	    y[jy] = y[jy] + temp1 * AP[kk];
	    ix = jx;
	    iy = jy;
	    for (k = kk + 1; k < kk + n - j; k++) {
		ix = ix + incx;
		iy = iy + incy;
		y[iy] = y[iy] + temp1 * AP[k];
		temp2 = temp2 + AP[k] * x[ix];
	    }
	    y[jy] = y[jy] + alpha * temp2;
	    jx = jx + incx;
	    jy = jy + incy;
	    kk = kk + n - j;
	}
    }
}
