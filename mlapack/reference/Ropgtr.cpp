/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ropgtr.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Ropgtr(const char *uplo, INTEGER n, REAL * ap, REAL * tau, REAL * q, INTEGER ldq, REAL * work, INTEGER * info)
{
    INTEGER i, j, ij;
    INTEGER iinfo;
    INTEGER upper;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -6;
    }
    if (*info != 0) {
	Mxerbla("Ropgtr", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
    if (upper) {
//Q was determined by a call to DSPTRD with UPLO = 'U'
//Unpack the vectors which define the elementary reflectors and
//set the last row and column of Q equal to those of the unit
//matrix
	ij = 2;
	for (j = 0; j < n - 1; j++) {
	    for (i = 0; i < j - 1; i++) {
		q[i + j * ldq] = ap[ij];
		ij++;

	    }
	    ij += 2;
	    q[n + j * ldq] = Zero;

	}
	for (i = 0; i < n - 1; i++) {
	    q[i + n * ldq] = Zero;

	}
	q[n + n * ldq] = One;
//Generate Q(1:n-1,1:n-1)
	Rorg2l(n - 1, n - 1, n - 1, &q[0], ldq, &tau[1], &work[0], &iinfo);
    } else {
//Q was determined by a call to DSPTRD with UPLO = 'L'.
//Unpack the vectors which define the elementary reflectors and
//set the first row and column of Q equal to those of the unit
//matrix
	q[ldq + 1] = One;
	for (i = 1; i < n; i++) {
	    q[i + ldq] = Zero;

	}
	ij = 3;
	for (j = 2; j <= n; j++) {
	    q[j * ldq + 1] = Zero;
	    for (i = j + 1; i < n; i++) {
		q[i + j * ldq] = ap[ij];
		ij++;
	    }
	    ij += 2;
	}
	if (n > 1) {
//Generate Q(2:n,2:n)
	    Rorg2r(n - 1, n - 1, n - 1, &q[(ldq * 2) + 2], ldq, &tau[0], &work[0], &iinfo);
	}
    }
    return;
}
