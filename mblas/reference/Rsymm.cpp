/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rsymm.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 *
 * $Id: Rsymm.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $

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

/*
Based on http://www.netlib.org/blas/dsymm.f
Rsymm performs one of the matrix-matrix operations
 C := alpha*A*B + beta*C,
 or
 C := alpha*B*A + beta*C,
where alpha and beta are scalars,  A is a symmetric matrix and  B and
C are m by n matrices.
*/

#include <mblas.h>

void Rsymm(const char *side, const char *uplo, INTEGER m, INTEGER n, REAL alpha, REAL * A, INTEGER lda, REAL * B, INTEGER ldb,
	   REAL beta, REAL * C, INTEGER ldc)
{
    INTEGER i, j, k, nrowa, upper, info;
    REAL Zero = 0.0, One = 1.0;
    REAL temp1, temp2;

//set nrowa as the number of rows of a.
    if (Mlsame(side, "L"))
	nrowa = m;
    else
	nrowa = n;
    upper = Mlsame(uplo, "U");

//test the input parameters.
    info = 0;
    if ((!Mlsame(side, "L")) && (!Mlsame(side, "R")))
	info = 1;
    else if ((!upper) && (!Mlsame(uplo, "L")))
	info = 2;
    else if (m < 0)
	info = 3;
    else if (n < 0)
	info = 4;
    else if (lda < max((INTEGER) 1, nrowa))
	info = 7;
    else if (ldb < max((INTEGER) 1, m))
	info = 9;
    else if (ldc < max((INTEGER) 1, m))
	info = 12;
    if (info != 0) {
	Mxerbla("Rsymm ", info);
	return;
    }
//quick return if possible.
    if ((m == 0) || (n == 0) || ((alpha == Zero) && (beta == One)))
	return;

//and when alpha==Zero.
    if (alpha == Zero) {
	if (beta == Zero) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = Zero;
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = beta * C[i + j * ldc];
		}
	    }
	}
	return;
    }
//start the operations.
    if (Mlsame(side, "L")) {
//Form C := alpha*A*B + beta*C.
	if (upper) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    temp1 = alpha * B[i + j * ldb];
		    temp2 = Zero;
		    for (k = 0; k < i; k++) {
			C[k + j * ldc] = C[k + j * ldc] + temp1 * A[k + i * lda];
			temp2 = temp2 + B[k + j * ldb] * A[k + i * lda];
		    }
		    if (beta == Zero) {
			C[i + j * ldc] = temp1 * A[i + i * lda] + alpha * temp2;
		    } else {
			C[i + j * ldc] = beta * C[i + j * ldc] + temp1 * A[i + i * lda] + alpha * temp2;
		    }
		}
	    }
	} else {
//Form  C := alpha*A*B + beta*C.
	    for (j = 0; j < n; j++) {
		for (i = m - 1; i >= 0; i--) {
		    temp1 = alpha * B[i + j * ldb];
		    temp2 = Zero;
		    for (k = i + 1; k < m; k++) {
			C[k + j * ldc] = C[k + j * ldc] + temp1 * A[k + i * lda];
			temp2 = temp2 + B[k + j * ldb] * A[k + i * lda];
		    }
		    if (beta == Zero) {
			C[i + j * ldc] = temp1 * A[i + i * lda] + alpha * temp2;
		    } else {
			C[i + j * ldc] = beta * C[i + j * ldc] + temp1 * A[i + i * lda] + alpha * temp2;
		    }
		}
	    }
	}
    } else {
//Form  C := alpha*B*A + beta*C.
	for (j = 0; j < n; j++) {
	    temp1 = alpha * A[j + j * lda];
	    if (beta == Zero) {
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = temp1 * B[i + j * ldb];
		}
	    } else {
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = beta * C[i + j * ldc] + temp1 * B[i + j * ldb];
		}
	    }
	    for (k = 0; k < j; k++) {
		if (upper)
		    temp1 = alpha * A[k + j * lda];
		else
		    temp1 = alpha * A[j + k * lda];

		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = C[i + j * ldc] + temp1 * B[i + k * ldb];
		}
	    }
	    for (k = j + 1; k < n; k++) {
		if (upper)
		    temp1 = alpha * A[j + k * lda];
		else
		    temp1 = alpha * A[k + j * lda];

		for (i = 0; i < m; i++) {
		    C[i + j * ldc] = C[i + j * ldc] + temp1 * B[i + k * ldb];
		}
	    }
	}
    }
    return;
}
