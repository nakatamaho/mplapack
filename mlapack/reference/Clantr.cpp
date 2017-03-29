/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clantr.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

REAL Clantr(const char *norm, const char *uplo, const char *diag, INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, REAL * work)
{
    INTEGER i, j;
    REAL sum, scale;
    INTEGER udiag;
    REAL value = 0.0;
    REAL mtemp1, mtemp2 = 0.0;
    REAL One = 1.0, Zero = 0.0;

    if (min(m, n) == 0) {
	value = Zero;
    } else if (Mlsame(norm, "M")) {

//Find max(abs(A(i,j))).
	if (Mlsame(diag, "U")) {
	    value = One;
	    if (Mlsame(uplo, "U")) {
		for (j = 0; j < n; j++) {
		    for (i = 0; i < min(m, j - 1); i++) {
			mtemp1 = value, mtemp2 = abs(A[i + j * lda]);
			value = max(mtemp1, mtemp2);
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (i = j + 1; i <= m; i++) {
			mtemp1 = value, mtemp1 = abs(A[i + j * lda]);
			value = max(mtemp1, mtemp2);
		    }
		}
	    }
	} else {
	    value = Zero;
	    if (Mlsame(uplo, "U")) {
		for (j = 0; j < n; j++) {
		    for (i = 0; i < min(m, j); i++) {
			mtemp1 = value, mtemp2 = abs(A[i + j * lda]);
			value = max(mtemp1, mtemp2);
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (i = j; i < m; i++) {
			mtemp1 = value, mtemp2 = abs(A[i + j * lda]);
			value = max(mtemp1, mtemp2);

		    }

		}
	    }
	}
    } else if (Mlsame(norm, "O") || Mlsame(norm, "1")) {
//Find norm1(A).

	value = Zero;
	udiag = Mlsame(diag, "U");
	if (Mlsame(uplo, "U")) {
	    for (j = 0; j < n; j++) {
		if (udiag && j <= m) {
		    sum = One;
		    for (i = 0; i < j - 1; i++) {
			sum += abs(A[i + j * lda]);
		    }
		} else {
		    sum = Zero;
		    for (i = 0; i < min(m, j); i++) {
			sum += abs(A[i + j * lda]);
		    }
		}
		value = max(value, sum);
	    }
	} else {
	    for (j = 0; j < n; j++) {
		if (udiag) {
		    sum = One;
		    for (i = j + 1; i < m; i++) {
			sum += abs(A[i + j * lda]);
		    }
		} else {
		    sum = Zero;
		    for (i = j; i < m; i++) {
			sum += abs(A[i + j * lda]);
		    }
		}
		value = max(value, sum);
	    }
	}
    } else if (Mlsame(norm, "I")) {
//Find normI(A).
	if (Mlsame(uplo, "U")) {
	    if (Mlsame(diag, "U")) {
		for (i = 0; i < m; i++) {
		    work[i] = One;
		}
		for (j = 0; j < n; j++) {
		    for (i = 0; i < min(m, j - 1); i++) {
			work[i] += abs(A[i + j * lda]);

		    }
		}
	    } else {
		for (i = 0; i < m; i++) {
		    work[i] = Zero;
		}
		for (j = 0; j < n; j++) {
		    for (i = 0; i < min(m, j); i++) {
			work[i] += abs(A[i + j * lda]);
		    }
		}
	    }
	} else {
	    if (Mlsame(diag, "U")) {
		for (i = 0; i < n; i++) {
		    work[i] = One;
		}
		for (i = n + 1; i < m; i++) {
		    work[i] = Zero;
		}
		for (j = 0; j < n; j++) {
		    for (i = j + 1; i < m; i++) {
			work[i] += abs(A[i + j * lda]);
		    }
		}
	    } else {
		for (i = 0; i < m; i++) {
		    work[i] = Zero;
		}
		for (j = 0; j < n; j++) {
		    for (i = j; i < m; i++) {
			work[i] += abs(A[i + j * lda]);
		    }
		}
	    }
	}
	value = Zero;
	for (i = 0; i < m; i++) {
	    mtemp1 = value, mtemp2 = work[i];
	    value = max(mtemp1, mtemp2);
	}
    } else if (Mlsame(norm, "F") || Mlsame(norm, "E")) {
//Find normF(A).
	if (Mlsame(uplo, "U")) {
	    if (Mlsame(diag, "U")) {
		scale = One;
		sum = min(m, n);
		for (j = 2; j < n; j++) {
		    Classq(min(m, j - 1), &A[j * lda], 1, &scale, &sum);

		}
	    } else {
		scale = Zero;
		sum = One;
		for (j = 0; j < n; j++) {
		    Classq(min(m, j), &A[j * lda], 1, &scale, &sum);

		}
	    }
	} else {
	    if (Mlsame(diag, "U")) {
		scale = One;
		sum = min(m, n);
		for (j = 0; j < n; j++) {
		    Classq(m - j, &A[min(m, j + 1) + j * lda], 1, &scale, &sum);
		}
	    } else {
		scale = Zero;
		sum = One;
		for (j = 0; j < n; j++) {
		    Classq(m - j + 1, &A[j + j * lda], 1, &scale, &sum);

		}
	    }
	}
	value = scale * sqrt(sum);
    }
    return value;
}
