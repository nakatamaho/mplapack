/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cher2.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Cher2.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/cher2.f
Cher2 performs the hermitian rank 2 operation
 A := alpha*x*conjg(y') + conjg(alpha)*y*conjg(x') + A,
where alpha is a scalar, x and y are n element vectors and A is an n
by n hermitian matrix.
*/

#include <mblas.h>

void Cher2(const char *uplo, INTEGER n, COMPLEX alpha, COMPLEX * x,
	   INTEGER incx, COMPLEX * y, INTEGER incy, COMPLEX * A, INTEGER lda)
{
    INTEGER i, ix, iy, j, jx, jy, kx, ky;
    REAL Zero = 0.0;
    COMPLEX temp1, temp2;

//test the input parameters.
    INTEGER info = 0;

    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L"))
	info = 1;
    else if (n < 0)
	info = 2;
    else if (incx == 0)
	info = 5;
    else if (incy == 0)
	info = 7;
    else if (lda < max((INTEGER) 1, n))
	info = 9;
    if (info != 0) {
	Mxerbla("Cher2 ", info);
	return;
    }
//quick return if possible.
    if ((n == 0) || (alpha == Zero))
	return;

    if (incx > 0)
	kx = 0;
    else
	kx = -(n - 1) * incx;
    if (incy > 0)
	ky = 0;
    else
	ky = -(n - 1) * incy;
    jx = kx;
    jy = ky;

    if (Mlsame(uplo, "U")) {
	for (j = 0; j < n; j++) {
	    if ((x[jx] != Zero) || (y[jy] != Zero)) {
		temp1 = alpha * conj(y[jy]);
		temp2 = conj(alpha * x[jx]);
		ix = kx;
		iy = ky;
		for (i = 0; i < j; i++) {
		    A[i + j * lda] = A[i + j * lda] + x[ix] * temp1 + y[iy] * temp2;
		    ix = ix + incx;
		    iy = iy + incy;
		}
		A[j + j * lda] = (A[j + j * lda]).real() + (x[jx] * temp1 + y[jy] * temp2).real();
	    } else {
		A[j + j * lda] = (A[j + j * lda]).real();
	    }
	    jx = jx + incx;
	    jy = jy + incy;
	}
    } else {
//form  a  when a is stored in the lower triangle.
	for (j = 0; j < n; j++) {
	    if ((x[jx] != Zero) || (y[jy] != Zero)) {
		temp1 = alpha * conj(y[jy]);
		temp2 = conj(alpha * x[jx]);
		A[j + j * lda] = (A[j + j * lda]).real() + (x[jx] * temp1 + y[jy] * temp2).real();
		ix = jx;
		iy = jy;
		for (i = j + 1; i < n; i++) {
		    ix = ix + incx;
		    iy = iy + incy;
		    A[i + j * lda] = A[i + j * lda] + x[ix] * temp1 + y[iy] * temp2;
		}
	    } else {
		A[j + j * lda] = (A[j + j * lda]).real();
	    }
	    jx = jx + incx;
	    jy = jy + incy;
	}
    }
    return;
}
