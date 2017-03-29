/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cptts2.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cptts2(INTEGER iuplo, INTEGER n, INTEGER nrhs, REAL * d, COMPLEX * e, COMPLEX * B, INTEGER ldb)
{
    INTEGER i, j;
    REAL One = 1.0;

    if (n <= 1) {
	if (n == 1) {
	    CRscal(nrhs, One / d[1], &B[0], ldb);
	}
	return;
    }
    if (iuplo == 1) {
//Solve A * X = B using the factorization A = U'*D*U,
//overwriting each right hand side vector with its solution.
	if (nrhs <= 2) {
	    j = 0;
	  L10:
//Solve U' * x = b.
	    for (i = 1; i < n; i++) {
		B[i + j * ldb] = B[i + j * ldb] - B[i - 1 + j * ldb] * conj(e[i - 1]);
	    }
//Solve D * U * x = b.
	    for (i = 0; i < n; i++) {
		B[i + j * ldb] = B[i + j * ldb] / d[i];
	    }
	    for (i = n - 1; i >= 1; i--) {
		B[i + j * ldb] = B[i + j * ldb] - B[i + 1 + j * ldb] * e[i];
	    }
	    if (j < nrhs) {
		j++;
		goto L10;
	    }
	} else {
	    for (j = 0; j < nrhs; j++) {
//Solve U' * x = b.
		for (i = 1; i < n; i++) {
		    B[i + j * ldb] = B[i + j * ldb] - B[i - 1 + j * ldb] * conj(e[i - 1]);
		}
//Solve D * U * x = b.
		B[n + j * ldb] = B[n + j * ldb] / d[n];
		for (i = n - 1; i >= 1; i--) {
		    B[i + j * ldb] = B[i + j * ldb] / d[i] - B[i + 1 + j * ldb] * e[i];
		}
	    }
	}
    } else {
//Solve A * X = B using the factorization A = L*D*L',
//overwriting each right hand side vector with its solution.
	if (nrhs <= 2) {
	    j = 0;
	  L80:
//Solve L * x = b.
	    for (i = 1; i < n; i++) {
		B[i + j * ldb] = B[i + j * ldb] - B[i - 1 + j * ldb] * e[i - 1];
	    }
//Solve D * L' * x = b.
	    for (i = 0; i < n; i++) {
		B[i + j * ldb] = B[i + j * ldb] / d[i];
	    }
	    for (i = n - 1; i >= 1; i--) {
		B[i + j * ldb] = B[i + j * ldb] - B[i + 1 + j * ldb] * conj(e[i]);
	    }
	    if (j < nrhs) {
		j++;
		goto L80;
	    }
	} else {
	    for (j = 0; j < nrhs; j++) {
//Solve L * x = b.
		for (i = 1; i < n; i++) {
		    B[i + j * ldb] = B[i + j * ldb] - B[i - 1 + j * ldb] * e[i - 1];
		}
//Solve D * L' * x = b.
		B[n + j * ldb] = B[n + j * ldb] / d[n];
		for (i = n - 1; i >= 1; i--) {
		    B[i + j * ldb] = B[i + j * ldb] / d[i] - B[i + 1 + j * ldb] * conj(e[i]);
		}
	    }
	}
    }
    return;
}
