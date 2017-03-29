/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clanhb.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

REAL Clanhb(const char *norm, const char *uplo, INTEGER n, INTEGER k, COMPLEX * AB, INTEGER ldab, REAL * work)
{
    INTEGER i, j, l;
    REAL sum, absa, scale;
    REAL value;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1;

    if (n == 0) {
	value = Zero;
    } else if (Mlsame(norm, "M")) {
//Find max(abs(A(i,j))).
	value = Zero;
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		for (i = max(k + 2 - j, (INTEGER) 1); i <= k; i++) {
		    value = max(value, abs(AB[i + j * ldab]));
		}
	    }
	    mtemp1 = abs(AB[k + 1 + j * ldab].real());
	    value = max(value, mtemp1);
	} else {
	    for (j = 0; j < n; j++) {
		mtemp1 = abs(AB[j * ldab + 1].real());
		value = max(value, mtemp1);
		for (i = 1; i < min(n + 1 - j, k + 1); i++) {
		    value = max(value, abs(AB[i + j * ldab]));
		}
	    }
	}
    } else if (Mlsame(norm, "I") || Mlsame(norm, "O") || Mlsame(norm, "1")) {
//Find normI(A) ( = norm1(A), since A is hermitian).
	value = Zero;
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		sum = Zero;
		l = k + 1 - j;
		for (i = max((INTEGER) 1, j - k); i <= j - 1; i++) {
		    absa = abs(AB[l + i + j * ldab]);
		    sum = sum + absa;
		    work[i] = work[i] + absa;
		}
		work[j] = sum + abs(AB[k + 1 + j * ldab].real());
	    }
	    for (i = 0; i < n; i++) {
		value = max(value, work[i]);
	    }
	} else {
	    for (i = 0; i < n; i++) {
		work[i] = Zero;
	    }
	    for (j = 0; j < n; j++) {
		sum = work[j] + abs(AB[j * ldab + 1].real());
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
		for (j = 2; j <= n; j++) {
		    Classq(min(j - 1, k), &AB[max(k + 2 - j, (INTEGER) 1) + j * ldab], 1, &scale, &sum);
		}
		l = k + 1;
	    } else {
		for (j = 0; j < n - 1; j++) {
		    Classq(min(n - j, k), &AB[j * ldab + 2], 1, &scale, &sum);
		}
		l = 0;
	    }
	    sum = sum * 2;
	} else {
	    l = 0;
	}
	Classq(n, &AB[l + ldab], ldab, &scale, &sum);
	value = scale * sqrt(sum);
	for (j = 0; j < n; j++) {
	    if (AB[l + j * ldab] != Zero) {
		absa = abs(AB[l + j * ldab].real());
		if (scale < absa) {
		    sum = sum * ((scale / absa) * (scale / absa)) + One;
		    scale = absa;
		} else {
		    sum = sum + (absa / scale) * (absa / scale);
		}
	    }

	}
    }
    return value;
}
