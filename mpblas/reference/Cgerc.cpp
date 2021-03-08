/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cgerc.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Cgerc.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/zgerc.f
Cgerc performs the rank 1 operation
 A := alpha*x*conjg( y' ) + A,
where alpha is a scalar, x is an m element vector, y is an n element
vector and A is an m by n matrix.
*/

#include <mpblas.h>

void Cgerc(INTEGER m, INTEGER n, COMPLEX alpha, COMPLEX * x, INTEGER incx, COMPLEX * y, INTEGER incy, COMPLEX * A, INTEGER lda)
{
    INTEGER i, ix, j, jy, kx;
    INTEGER info = 0;
    REAL Zero = 0.0;
    COMPLEX temp;

//Test the input parameters.
    if (m < 0)
	info = 1;
    else if (n < 0)
	info = 2;
    else if (incx == 0)
	info = 5;
    else if (incy == 0)
	info = 7;
    else if (lda < max((INTEGER) 1, m))
	info = 9;
    if (info != 0) {
	Mxerbla("Cgerc ", info);
	return;
    }
//quick return if possible.
    if ((m == 0) || (n == 0) || (alpha == Zero))
	return;

//start the operations. in this version the elements of a are
//accessed sequentially with one pass through A.
    if (incy > 0)
	jy = 0;
    else
	jy = (1 - n) * incy;

    if (incx > 0)
	kx = 0;
    else
	kx = (1 - m) * incx;

    for (j = 0; j < n; j++) {
	if (y[jy] != Zero) {
	    temp = alpha * conj(y[jy]);
	    ix = kx;
	    for (i = 0; i < m; i++) {
		A[i + j * lda] = A[i + j * lda] + x[ix] * temp;
		ix = ix + incx;
	    }
	}
	jy = jy + incy;
    }
    return;
}
