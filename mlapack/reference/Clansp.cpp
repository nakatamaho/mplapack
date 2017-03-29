/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clansp.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

REAL Clansp(const char *norm, const char *uplo, INTEGER n, COMPLEX * ap, REAL * work)
{
    INTEGER i, j, k;
    REAL sum, absa, scale;
    REAL value = 0.0;
    REAL mtemp1, mtemp2;
    REAL Zero = 0.0, One = 1.0;

    if (n == 0) {
	value = Zero;
    } else if (Mlsame(norm, "M")) {
//Find max(abs(A(i,j))).
	value = Zero;
	if (Mlsame(uplo, "U")) {
	    k = 0;
	    for (j = 0; j < n; j++) {
		for (i = k; i <= k + j - 1; i++) {
		    mtemp1 = value, mtemp2 = abs(ap[i]);
		    value = max(mtemp1, mtemp2);
		}
		k += j;
	    }
	} else {
	    k = 0;
	    for (j = 0; j < n; j++) {
		for (i = k; i <= k + n - j; i++) {
		    mtemp1 = value, mtemp2 = abs(ap[i]);
		    value = max(mtemp1, mtemp2);
		}
		k = k + n - j + 1;
	    }
	}
    } else if (Mlsame(norm, "I") || Mlsame(norm, "O") || Mlsame(norm, "1")) {
//Find normI(A) ( = norm1(A), since A is symmetric).
	value = Zero;
	k = 0;
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		sum = Zero;
		for (i = 0; i < j - 1; i++) {
		    absa = abs(ap[k]);
		    sum = sum + absa;
		    work[i] = work[i] + absa;
		    k++;
		}
		work[j] = sum + abs(ap[k]);
		k++;
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
		sum = work[j] + abs(ap[k]);
		k++;
		for (i = j + 1; i <= n; i++) {
		    absa = abs(ap[k]);
		    sum = sum + absa;
		    work[i] = work[i] + absa;
		    k++;
		}
		value = max(value, sum);
	    }
	}
    } else if (Mlsame(norm, "F") || Mlsame(norm, "E")) {
//Find normF(A).
	scale = Zero;
	sum = One;
	k = 2;
	if (Mlsame(uplo, "U")) {
	    for (j = 2; j <= n; j++) {
		Classq(j - 1, &ap[k], 1, &scale, &sum);
		k += j;
	    }
	} else {
	    for (j = 0; j < n - 1; j++) {
		Classq(n - j, &ap[k], 1, &scale, &sum);
		k = k + n - j + 1;

	    }
	}
	sum = sum * 2;
	k = 0;
	for (i = 0; i < n; i++) {
	    if (ap[k].real() != Zero) {
		absa = abs(ap[k].real());
		if (scale < absa) {
		    sum = sum * ((scale / absa) * (scale / absa)) + One;
		    scale = absa;
		} else {
		    sum = sum + (absa / scale) * (absa / scale);
		}
	    }
	    if (ap[k].imag() != Zero) {
		absa = abs(ap[k]);
		if (scale < absa) {
		    sum = sum * ((scale / absa) * (scale / absa)) + One;
		    scale = absa;
		} else {
		    sum = sum + (absa / scale) * (absa / scale);
		}
	    }
	    if (Mlsame(uplo, "U")) {
		k = k + i + 1;
	    } else {
		k = k + n - i + 1;
	    }

	}
	value = scale * sqrt(sum);
    }
    return value;
}
