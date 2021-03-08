/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemv.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Rgemv.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/dgemv.f
Rgemv performs one of the matrix-vector operations
y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,
where alpha and beta are scalars, x and y are vectors and A is an
m by n matrix.
*/

#include <mpblas.h>

void Rgemv(const char *trans, INTEGER m, INTEGER n, REAL alpha, REAL * A, INTEGER lda, REAL * x, INTEGER incx, REAL beta, REAL * y,
	   INTEGER incy)
{
    INTEGER lenx, leny, i, ix, jx, kx, iy, j, jy, ky;
    INTEGER info = 0;
    REAL Zero = 0.0, One = 1.0;
    REAL temp;

//Test the input parameters.
    if (!Mlsame(trans, "N") && !Mlsame(trans, "T") && !Mlsame(trans, "C"))
	info = 1;
    else if (m < 0)
	info = 2;
    else if (n < 0)
	info = 3;
    else if (lda < max((INTEGER) 1, m))
	info = 6;
    else if (incx == 0)
	info = 8;
    else if (incy == 0)
	info = 11;
    if (info != 0) {
	Mxerbla("Rgemv ", info);
	return;
    }
//Quick return if possible.
    if ((m == 0) || (n == 0) || ((alpha == Zero) && (beta == One)))
	return;

//Set lenx and leny, the lengths of the vectors x and y, and set
//up the start points in x and y.
    if (Mlsame(trans, "N")) {
	lenx = n;
	leny = m;
    } else {
	lenx = m;
	leny = n;
    }
    if (incx > 0)
	kx = 0;
    else
	kx = (1 - lenx) * incx;
    if (incy > 0)
	ky = 0;
    else
	ky = (1 - leny) * incy;

//start the operations. in this version the elements of a are
//accessed sequentially with One pass through a.
//first form  y := beta*y.
    if (beta != One) {
	iy = ky;
	if (beta == Zero) {
	    for (i = 0; i < leny; i++) {
		y[iy] = Zero;
		iy = iy + incy;
	    }
	} else {
	    for (i = 0; i < leny; i++) {
		y[iy] = beta * y[iy];
		iy = iy + incy;
	    }
	}
    }
    if (alpha == Zero)
	return;
    if (Mlsame(trans, "N")) {
//form y := alpha*A*x + y.
	jx = kx;
	for (j = 0; j < n; j++) {
	    if (x[jx] != Zero) {
		temp = alpha * x[jx];
		iy = ky;
		for (i = 0; i < m; i++) {
		    y[iy] = y[iy] + temp * A[i + j * lda];
		    iy = iy + incy;
		}
	    }
	    jx = jx + incx;
	}
    } else {
//Form y := alpha*A'*x + y.
	jy = ky;
	for (j = 0; j < n; j++) {
	    temp = Zero;
	    ix = kx;
	    for (i = 0; i < m; i++) {
		temp = temp + A[i + j * lda] * x[ix];
		ix = ix + incx;
	    }
	    y[jy] = y[jy] + alpha * temp;
	    jy = jy + incy;
	}
    }
    return;
}
