/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlange.cpp,v 1.6 2010/08/12 22:52:45 nakatamaho Exp $ 
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

#include <mpblas.h>
#include <mplapack.h>

REAL Rlange(const char *norm, INTEGER m, INTEGER n, REAL * A, INTEGER lda, REAL * work)
{
    INTEGER i, j;
    REAL sum, scale;
    REAL value = 0.0;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp, mtemp1;

    if (min(m, n) == 0) {
	value = Zero;
    } else if (Mlsame(norm, "M")) {

//Find max(abs(A(i,j))).
	value = Zero;
	for (j = 0; j < n; j++) {
	    for (i = 0; i < m; i++) {
		mtemp = value;
		mtemp1 = abs(A[i + j * lda]);
		value = max(mtemp, mtemp1);
	    }
	}
    } else if (Mlsame(norm, "O") || (Mlsame(norm, "1"))) {
//Find norm1(A).
	value = Zero;
	for (j = 0; j < n; j++) {
	    sum = Zero;
	    for (i = 0; i < m; i++) {
		sum += abs(A[i + j * lda]);
	    }
            mtemp = value;
	    mtemp1 = sum;
	    value = max(mtemp, mtemp1);
	}
    } else if (Mlsame(norm, "I")) {
//Find normI(A).
	for (i = 0; i < m; i++) {
	    work[i] = Zero;
	}
	for (j = 0; j < n; j++) {
	    for (i = 0; i < m; i++) {
		work[i] += abs(A[i + j * lda]);
	    }
	}
	value = Zero;
	for (i = 0; i < m; i++) {
	    value = max(value, work[i]);
	}
    } else if (Mlsame(norm, "F") || Mlsame(norm, "E")) {
//Find normF(A).
	scale = Zero;
	sum = One;
	for (j = 0; j < n; j++) {
	    Rlassq(m, &A[j * lda], 1, &scale, &sum);
	}
	value = scale * sqrt(sum);
    }
    return value;
}
