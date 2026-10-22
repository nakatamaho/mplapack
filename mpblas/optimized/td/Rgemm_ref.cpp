/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemm.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Rgemm.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/dgemm.f
Rgemm performs one of the matrix-matrix operations
 C := alpha*op(A)*op(B) + beta*C,
where op(X) is one of
 op(X) = X or op(X) = X',
alpha and beta are scalars, and A, B and C are matrices, with op( A )
an m by k matrix, op(B) a k by n matrix and C an m by n matrix.
*/

#include <mpblas.h>

void Rgemm_ref(const char *transa, const char *transb, INTEGER m, INTEGER n, INTEGER k, REAL alpha, REAL * A, INTEGER lda, REAL * B,
	   INTEGER ldb, REAL beta, REAL * C, INTEGER ldc)
{
    INTEGER i, j, l, nota, notb, nrowa, ncola, nrowb, info;
    REAL temp;
    REAL Zero = 0.0, One = 1.0;

    nota = Mlsame(transa, "N");
    notb = Mlsame(transb, "N");
    if (nota) {
	nrowa = m;
	ncola = k;
    } else {
	nrowa = k;
	ncola = m;
    }
    if (notb) {
	nrowb = k;
    } else {
	nrowb = n;
    }
//Test the input parameters.
    info = 0;
    if (!nota && (!Mlsame(transa, "C")) && (!Mlsame(transa, "T")))
	info = 1;
    else if (!notb && (!Mlsame(transb, "C")) && (!Mlsame(transb, "T")))
	info = 2;
    else if (m < 0)
	info = 3;
    else if (n < 0)
	info = 4;
    else if (k < 0)
	info = 5;
    else if (lda < max((INTEGER) 1, nrowa))
	info = 8;
    else if (ldb < max((INTEGER) 1, nrowb))
	info = 10;
    else if (ldc < max((INTEGER) 1, m))
	info = 13;
    if (info != 0) {
	Mxerbla("Rgemm ", info);
	return;
    }
//Quick return if possible.
    if ((m == 0) || (n == 0) || (((alpha == Zero) || (k == 0)) && (beta == One)))
	return;

//And when alpha == 0.0
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
//Start the operations.
    if (notb) {
	if (nota) {
//Form C := alpha*A*B + beta*C.
	    for (j = 0; j < n; j++) {
		if (beta == Zero) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] = Zero;
		    }
		} else if (beta != One) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
		for (l = 0; l < k; l++) {
		    if (B[l + j * ldb] != Zero) {
			temp = alpha * B[l + j * ldb];
			for (i = 0; i < m; i++) {
			    C[i + j * ldc] = C[i + j * ldc] + temp * A[i + l * lda];
			}
		    }
		}
	    }
	} else {
//Form  C := alpha*A'*B + beta*C.
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    temp = Zero;
		    for (l = 0; l < k; l++) {
			temp = temp + A[l + i * lda] * B[l + j * ldb];
		    }
		    if (beta == Zero)
			C[i + j * ldc] = alpha * temp;
		    else
			C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
		}
	    }
	}
    } else {
	if (nota) {
//Form  C := alpha*A*B' + beta*C.
	    for (j = 0; j < n; j++) {
		if (beta == Zero) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] = Zero;
		    }
		} else if (beta != One) {
		    for (i = 0; i < m; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
		for (l = 0; l < k; l++) {
		    if (B[j + l * ldb] != Zero) {
			temp = alpha * B[j + l * ldb];
			for (i = 0; i < m; i++) {
			    C[i + j * ldc] = C[i + j * ldc] + temp * A[i + l * lda];
			}
		    }
		}
	    }
	} else {
//Form  C := alpha*A'*B' + beta*C.
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    temp = Zero;
		    for (l = 0; l < k; l++) {
			temp = temp + A[l + i * lda] * B[j + l * ldb];
		    }
		    if (beta == Zero)
			C[i + j * ldc] = alpha * temp;
		    else
			C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
		}
	    }
	}
    }
    return;
}
