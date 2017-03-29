/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Chpr.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Chpr.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/chpr.f
Chpr performs the hermitian rank 1 operation
 A := alpha*x*conjg(x') + A,
where alpha is a real scalar, x is an n element vector and A is an
n by n hermitian matrix, supplied in packed form.
*/

#include <mblas.h>

void Chpr(const char *uplo, INTEGER n, REAL alpha, COMPLEX * x, INTEGER incx, COMPLEX * AP)
{
    INTEGER ix, j, jx, k, kx, kk;
    INTEGER info = 0;
    REAL Zero = 0.0;
    COMPLEX temp;

//test the input parameters.
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L"))
	info = 1;
    else if (n < 0)
	info = 2;
    else if (incx == 0)
	info = 5;
    if (info != 0) {
	Mxerbla("Chpr  ", info);
	return;
    }
//quick return if possible.
    if ((n == 0) || (alpha == Zero))
	return;

//set the start point in x if the increment is not unity.
    if (incx <= 0)
	kx = -(n - 1) * incx;
    else
	kx = 0;

//start the operations. in this version the elements of the array ap
//are accessed sequentially with one pass through ap.
    kk = 0;
    if (Mlsame(uplo, "U")) {
//Form  A  when upper triangle is stored in AP.
	jx = kx;
	for (j = 0; j < n; j++) {
	    if (x[jx] != Zero) {
		temp = alpha * conj(x[jx]);
		ix = kx;
		for (k = kk; k < kk + j; k++) {
		    AP[k] = AP[k] + x[ix] * temp;
		    ix = ix + incx;
		}
		AP[kk + j] = (AP[kk + j]).real() + (x[jx] * temp).real();
	    } else {
		AP[kk + j] = (AP[kk + j]).real();
	    }
	    jx = jx + incx;
	    kk = kk + j + 1;
	}
    } else {
//Form  A  when lower triangle is stored in AP.
	jx = kx;
	for (j = 0; j < n; j++) {
	    if (x[jx] != Zero) {
		temp = alpha * conj(x[jx]);
		AP[kk] = (AP[kk]).real() + (temp * x[jx]).real();
		ix = jx;
		for (k = kk + 1; k < kk + n - j; k++) {
		    ix = ix + incx;
		    AP[k] = AP[k] + x[ix] * temp;
		}
	    } else {
		AP[kk] = (AP[kk]).real();
	    }
	    jx = jx + incx;
	    kk = kk + n - j;
	}
    }
}
