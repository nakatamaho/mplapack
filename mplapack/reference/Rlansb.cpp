/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlansb.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

REAL Rlansb(const char *norm, const char *uplo, INTEGER n, INTEGER k, REAL * AB, INTEGER ldab, REAL * work)
{
    INTEGER i, j, l;
    REAL sum, absa, scale;
    REAL value = 0.0;
    REAL mtemp1, mtemp2 = 0.0;
    REAL One = 1.0, Zero = 0.0;

    if (n == 0) {
	value = Zero;
    } else if (Mlsame(norm, "M")) {
//Find max(abs(A(i,j))).
	value = Zero;
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		for (i = max(k + 2 - j, (INTEGER) 1); i <= k + 1; i++) {
		    mtemp1 = value, mtemp2 = abs(AB[i + j * ldab]);
		    value = max(mtemp1, mtemp2);
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < min(n + 1 - j, k + 1); i++) {
		    mtemp1 = value, mtemp1 = abs(AB[i + j * ldab]);
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
		l = k + 1 - j;
		for (i = max((INTEGER) 1, j - k); i <= j - 1; i++) {
		    absa = abs(AB[l + i + j * ldab]);
		    sum += absa;
		    work[i] += absa;
		}
		work[j] = sum + abs(AB[k + 1 + j * ldab]);
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
		sum = work[j] + abs(AB[j * ldab + 1]);
		l = 0 - j;
		for (i = j + 1; i <= min(n, j + k); i++) {
		    absa = abs(AB[l + i + j * ldab]);
		    sum += absa;
		    work[i] += absa;
		}
		value = max(value, sum);
	    }
	}
    } else if (Mlsame(norm, "F") || Mlsame(norm, "E")) {
//Find normF(A).
	scale = Zero;
	sum = One;
	if (k > 0) {
	    if (Mlsame(uplo, "U")) {
		for (j = 2; j < n; j++) {
		    Rlassq(min(j - 1, k), &AB[max(k + 2 - j, (INTEGER) 1) + j * ldab], 1, &scale, &sum);
		}
		l = k + 1;
	    } else {
		for (j = 0; j < n - 1; j++) {
		    Rlassq(min(n - j, k), &AB[j * ldab + 2], 1, &scale, &sum);

		}
		l = 0;
	    }
	    sum *= 2;
	} else {
	    l = 0;
	}
	Rlassq(n, &AB[l + ldab], ldab, &scale, &sum);
	value = scale * sqrt(sum);
    }
    return value;
}
