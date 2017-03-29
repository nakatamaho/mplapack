/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgtsv.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgtsv(INTEGER n, INTEGER nrhs, REAL * dl, REAL * d, REAL * du, REAL * B, INTEGER ldb, INTEGER * info)
{
    INTEGER i, j;
    REAL fact, temp;
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
	Mxerbla("Rgtsv ", -(*info));
	return;
    }

    if (n == 0) {
	return;
    }

    if (nrhs == 1) {

	for (i = 0; i < n - 2; i++) {
	    if (abs(d[i]) >= abs(dl[i])) {
//No row interchange required
		if (d[i] != Zero) {
		    fact = dl[i] / d[i];
		    d[i + 1] -= fact * du[i];
		    B[i + 1 + ldb] -= fact * B[i + ldb];
		} else {
		    *info = i;
		    return;
		}
		dl[i] = Zero;
	    } else {
//Interchange rows I and I+1
		fact = d[i] / dl[i];
		d[i] = dl[i];
		temp = d[i + 1];
		d[i + 1] = du[i] - fact * temp;
		dl[i] = du[i + 1];
		du[i + 1] = -fact * dl[i];
		du[i] = temp;
		temp = B[i + ldb];
		B[i + ldb] = B[i + 1 + ldb];
		B[i + 1 + ldb] = temp - fact * B[i + 1 + ldb];
	    }

	}
	if (n > 1) {
	    i = n - 1;
	    if (abs(d[i]) >= abs(dl[i])) {
		if (d[i] != Zero) {
		    fact = dl[i] / d[i];
		    d[i + 1] -= fact * du[i];
		    B[i + 1 + ldb] -= fact * B[i + ldb];
		} else {
		    *info = i;
		    return;
		}
	    } else {
		fact = d[i] / dl[i];
		d[i] = dl[i];
		temp = d[i + 1];
		d[i + 1] = du[i] - fact * temp;
		du[i] = temp;
		temp = B[i + ldb];
		B[i + ldb] = B[i + 1 + ldb];
		B[i + 1 + ldb] = temp - fact * B[i + 1 + ldb];
	    }
	}
	if (d[n] == Zero) {
	    *info = n;
	    return;
	}
    } else {
	for (i = 0; i < n - 2; i++) {
	    if (abs(d[i]) >= abs(dl[i])) {
//No row interchange required

		if (d[i] != Zero) {
		    fact = dl[i] / d[i];
		    d[i + 1] -= fact * du[i];
		    for (j = 0; j < nrhs; j++) {
			B[i + 1 + j * ldb] -= fact * B[i + j * ldb];

		    }
		} else {
		    *info = i;
		    return;
		}
		dl[i] = Zero;
	    } else {

//Interchange rows I and I+1
		fact = d[i] / dl[i];
		d[i] = dl[i];
		temp = d[i + 1];
		d[i + 1] = du[i] - fact * temp;
		dl[i] = du[i + 1];
		du[i + 1] = -fact * dl[i];
		du[i] = temp;
		for (j = 0; j < nrhs; j++) {
		    temp = B[i + j * ldb];
		    B[i + j * ldb] = B[i + 1 + j * ldb];
		    B[i + 1 + j * ldb] = temp - fact * B[i + 1 + j * ldb];

		}
	    }

	}
	if (n > 1) {
	    i = n - 1;
	    if (abs(d[i]) >= abs(dl[i])) {
		if (d[i] != Zero) {
		    fact = dl[i] / d[i];
		    d[i + 1] -= fact * du[i];
		    for (j = 0; j < nrhs; j++) {
			B[i + 1 + j * ldb] -= fact * B[i + j * ldb];

		    }
		} else {
		    *info = i;
		    return;
		}
	    } else {
		fact = d[i] / dl[i];
		d[i] = dl[i];
		temp = d[i + 1];
		d[i + 1] = du[i] - fact * temp;
		du[i] = temp;

		for (j = 0; j < nrhs; j++) {
		    temp = B[i + j * ldb];
		    B[i + j * ldb] = B[i + 1 + j * ldb];
		    B[i + 1 + j * ldb] = temp - fact * B[i + 1 + j * ldb];

		}
	    }
	}
	if (d[n] == Zero) {
	    *info = n;
	    return;
	}
    }

//Back solve with the matrix U from the factorization.
    if (nrhs <= 2) {
	j = 0;
	while (1) {
	    B[n + j * ldb] /= d[n];
	    if (n > 1) {
		B[n - 1 + j * ldb] = (B[n - 1 + j * ldb] - du[n - 1] * B[n + j * ldb]) / d[n - 1];
	    }
	    for (i = n - 2; i >= 1; i--) {
		B[i + j * ldb] = (B[i + j * ldb] - du[i] * B[i + 1 + j * ldb] - dl[i] * B[i + 2 + j * ldb]) / d[i];

	    }
	    if (j >= nrhs)
		break;
	    j++;
	    continue;
	}
    } else {

	for (j = 0; j < nrhs; j++) {
	    B[n + j * ldb] /= d[n];
	    if (n > 1) {
		B[n - 1 + j * ldb] = (B[n - 1 + j * ldb] - du[n - 1]
				      * B[n + j * ldb]) / d[n - 1];
	    }
	    for (i = n - 2; i >= 1; i--) {
		B[i + j * ldb] = (B[i + j * ldb] - du[i] * B[i + 1 + j * ldb] - dl[i] * B[i + 2 + j * ldb])
		    / d[i];

	    }

	}
    }
    return;
}
