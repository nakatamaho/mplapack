/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlansy.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

REAL Rlansy(const char *norm, const char *uplo, INTEGER n, REAL * A, INTEGER lda, REAL * work)
{
    INTEGER i, j;
    REAL sum, absa, scale, value = 0.0;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

    if (n == 0) {
	value = Zero;
	return value;
    }
    if (Mlsame(norm, "M")) {
//Find max(abs(A(i,j))).
	value = Zero;
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= j; i++) {
		    mtemp1 = value, mtemp2 = abs(A[i + j * lda]);
		    value = max(mtemp1, mtemp2);
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
		    mtemp1 = value, mtemp2 = abs(A[i + j * lda]);
		    value = max(mtemp1, mtemp2);
		}
	    }
	}
    } else if (Mlsame(norm, "I") || Mlsame(norm, "O") || Mlsame(norm, "1")) {
//Find normI(A) ( = norm1(A), since A is symmetric).
	value = Zero;
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		sum = Zero;
		for (i = 0; i < j; i++) {
		    absa = abs(A[i + j * lda]);
		    sum = sum + absa;
		    work[i] = work[i] + absa;
		}
		work[j] = sum + abs(A[j + j * lda]);
	    }
	    for (i = 0; i < n; i++) {
		mtemp1 = value, mtemp2 = work[i];
		value = max(mtemp1, mtemp2);
	    }
	} else {
	    for (i = 0; i < n; i++) {
		work[i] = Zero;
	    }
	    for (j = 0; j < n; j++) {
		sum = work[j] + abs(A[j + j * lda]);
		for (i = j + 1; i < n; i++) {
		    absa = abs(A[i + j * lda]);
		    sum = sum + absa;
		    work[i] = work[i] + absa;
		}
		mtemp1 = value, mtemp2 = sum;
		value = max(mtemp1, mtemp2);
	    }
	}
    } else if (Mlsame(norm, "F") || Mlsame(norm, "E")) {
//Find normF(A).
	scale = Zero;
	sum = One;
	if (Mlsame(uplo, "U")) {
	    for (j = 1; j < n; j++) {
		Rlassq(j, &A[j * lda], 1, &scale, &sum);
	    }
	} else {
	    for (j = 0; j < n - 1; j++) {
		Rlassq(n - j - 1, &A[(j + 1) + j * lda], 1, &scale, &sum);
	    }
	}
	sum = sum * Two;
	Rlassq(n, A, lda + 1, &scale, &sum);
	value = scale * sqrt(sum);
    }
    return value;
}
