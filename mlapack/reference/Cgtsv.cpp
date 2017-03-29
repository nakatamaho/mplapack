/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgtsv.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>

void Cgtsv(INTEGER n, INTEGER nrhs, COMPLEX * dl, COMPLEX * d, COMPLEX * du, COMPLEX * B, INTEGER ldb, INTEGER * info)
{
    INTEGER j, k;
    COMPLEX temp, mult;
    REAL Zero = 0.0;

    *info = 0;
    if (n < 0) {
	*info = -1;
    } else if (nrhs < 0) {
	*info = -2;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -7;
    }
    if (*info != 0) {
	Mxerbla("Cgtsv ", -(*info));
	return;
    }
    if (n == 0) {
	return;
    }
    for (k = 0; k < n - 1; k++) {
	if (dl[k] == Zero) {
//Subdiagonal is zero, no elimination is required.
	    if (d[k] == Zero) {
//Diagonal is zero: set INFO = K and return; a unique
//solution can not be found.
		*info = k;
		return;
	    }
	} else if (abs(d[k].real()) + abs(d[k].imag()) >= abs(dl[k].real()) + abs(dl[k].imag())) {
//No row INTEGERerchange required
	    mult = dl[k] / d[k];
	    d[k + 1] = d[k + 1] - mult * du[k];
	    for (j = 0; j < nrhs; j++) {
		B[k + 1 + j * ldb] = B[k + 1 + j * ldb] - mult * B[k + j * ldb];
	    }
	    if (k < n - 1) {
		dl[k] = Zero;
	    }
	} else {
//interchange rows K and K+1
	    mult = d[k] / dl[k];
	    d[k] = dl[k];
	    temp = d[k + 1];
	    d[k + 1] = du[k] - mult * temp;
	    if (k < n - 1) {
		dl[k] = du[k + 1];
		du[k + 1] = -mult * dl[k];
	    }
	    du[k] = temp;
	    for (j = 0; j < nrhs; j++) {
		temp = B[k + j * ldb];
		B[k + j * ldb] = B[k + 1 + j * ldb];
		B[k + 1 + j * ldb] = temp - mult * B[k + 1 + j * ldb];
	    }
	}
    }
    if (d[n] == Zero) {
	*info = n;
	return;
    }
//Back solve with the matrix U from the factorization.
    for (j = 0; j < nrhs; j++) {
	B[n + j * ldb] = B[n + j * ldb] / d[n];
	if (n > 1) {
	    B[n - 1 + j * ldb] = (B[n - 1 + j * ldb] - du[n - 1] * B[n + j * ldb]) / d[n - 1];
	}
	for (k = n - 2; k >= 1; k--) {
	    B[k + j * ldb] = (B[k + j * ldb] - du[k] * B[k + 1 + j * ldb] - dl[k] * B[k + 2 + j * ldb]) / d[k];

	}
    }
    return;
}
