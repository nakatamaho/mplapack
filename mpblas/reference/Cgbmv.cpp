/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cgbmv.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Cgbmv.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/zgbmv.f
Cgbmv performs one of the matrix-vector operations
y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, or
y := alpha*conjg(A')*x + beta*y
*/

#include <mpblas.h>

void Cgbmv(const char *trans, INTEGER m, INTEGER n, INTEGER kl, INTEGER ku,
	   COMPLEX alpha, COMPLEX * A, INTEGER lda, COMPLEX * x, INTEGER incx, COMPLEX beta, COMPLEX * y, INTEGER incy)
{
    INTEGER lenx, leny, i, ix, iy, j, jx, jy, k, kx, ky, noconj;
    REAL Zero = 0.0, One = 1.0;
    COMPLEX temp;

//Test the input parameters.
    int info = 0;
    if (!Mlsame(trans, "N") && !Mlsame(trans, "T") && !Mlsame(trans, "C"))
	info = 1;
    else if (m < 0)
	info = 2;
    else if (n < 0)
	info = 3;
    else if (kl < 0)
	info = 4;
    else if (ku < 0)
	info = 5;
    else if (lda < (kl + ku + 1))
	info = 8;
    else if (incx == 0)
	info = 10;
    else if (incy == 0)
	info = 13;

    if (info != 0) {
	Mxerbla("Cgbmv ", info);
	return;
    }
    //quick return if possible.
    if ((m == 0) || (n == 0) || (alpha == Zero && beta == One))
	return;

    noconj = Mlsame(trans, "T");

    //Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    //up the start points in  X  and  Y.
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

//start the operations. in this version the elements of A are
//accessed sequentially with one pass through the band part of A.
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
//Form y := alpha*A*x + y.
	jx = kx;
	for (j = 0; j < n; j++) {
	    if (x[jx] != Zero) {
		temp = alpha * x[jx];
		iy = ky;
		k = ku - j;
		for (i = max((INTEGER) 0, j - ku); i < min(m, j + kl + 1); i++) {
		    y[iy] = y[iy] + temp * A[(k + i) + j * lda];
		    iy = iy + incy;
		}
	    }
	    jx = jx + incx;
	    if (j >= ku)
		ky = ky + incy;
	}
    } else {
//Form y:=alpha*A'*x+y  or  y:=alpha*conjg(A')*x+y.
	jy = ky;
	for (j = 0; j < n; j++) {
	    temp = Zero;
	    ix = kx;
	    k = ku - j;
	    if (noconj) {
		for (i = max((INTEGER) 0, j - ku); i < min(m, j + kl + 1); i++) {
		    temp = temp + A[k + i + j * lda] * x[ix];
		    ix = ix + incx;
		}
	    } else {
		for (i = max((INTEGER) 0, j - ku); i < min(m, j + kl + 1); i++) {
		    temp = temp + conj(A[k + i + j * lda]) * x[ix];
		    ix = ix + incx;
		}
	    }
	    y[jy] = y[jy] + alpha * temp;
	    jy = jy + incy;
	    if (j >= ku)
		kx = kx + incx;
	}
    }
}
