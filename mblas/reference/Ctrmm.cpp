/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Ctrmm.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Ctrmm.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $

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
Based on http://www.netlib.org/blas/ctrmm.f
Ctrmm performs one of the matrix-matrix operations
 B := alpha*op(A)*B, or B := alpha*B*op(A)
where alpha is a scalar, B is an m by n matrix, A is a unit, or
non-unit, upper or lower triangular matrix and op(A) is one  of
op(A)=A or op(A)=A' or op(A)=conjg(A').
*/

#include <mblas.h>

void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER m, INTEGER n, COMPLEX alpha,
	   COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb)
{
    INTEGER i, info, j, k, lside, nrowa, nounit, upper;
    INTEGER noconj;
    REAL One = 1.0, Zero = 0.0;
    COMPLEX temp;

//test the input parameters.
    lside = Mlsame(side, "L");
    if (lside)
	nrowa = m;
    else
	nrowa = n;

    noconj = Mlsame(transa, "T");
    nounit = Mlsame(diag, "N");
    upper = Mlsame(uplo, "U");
    info = 0;
    if ((!lside) && (!Mlsame(side, "R")))
	info = 1;
    else if ((!upper) && (!Mlsame(uplo, "L")))
	info = 2;
    else if ((!Mlsame(transa, "N")) && (!Mlsame(transa, "T")) && (!Mlsame(transa, "C")))
	info = 3;
    else if ((!Mlsame(diag, "U")) && (!Mlsame(diag, "N")))
	info = 4;
    else if (m < 0)
	info = 5;
    else if (n < 0)
	info = 6;
    else if (lda < max((INTEGER) 1, nrowa))
	info = 9;
    else if (ldb < max((INTEGER) 1, m))
	info = 11;
    if (info != 0) {
	Mxerbla("Ctrmm ", info);
	return;
    }
//quick return if possible.
    if (m == 0 || n == 0)
	return;

//and when alpha==Zero.
    if (alpha == Zero) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < m; i++) {
		B[i + j * ldb] = Zero;
	    }
	}
	return;
    }
//start the operations.
    if (lside) {
	if (Mlsame(transa, "N")) {
//Form B := alpha*A*B.
	    if (upper) {
		for (j = 0; j < n; j++) {
		    for (k = 0; k < m; k++) {
			if (B[k + j * ldb] != Zero) {
			    temp = alpha * B[k + j * ldb];
			    for (i = 0; i < k; i++) {
				B[i + j * ldb] = B[i + j * ldb] + temp * A[i + k * lda];
			    }
			    if (nounit)
				temp = temp * A[k + k * lda];
			    B[k + j * ldb] = temp;
			}
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (k = m - 1; k >= 0; k--) {
			if (B[k + j * ldb] != Zero) {
			    temp = alpha * B[k + j * ldb];
			    B[k + j * ldb] = temp;
			    if (nounit)
				B[k + j * ldb] = B[k + j * ldb] * A[k + k * lda];
			    for (i = k + 1; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] + temp * A[i + k * lda];
			    }
			}
		    }
		}
	    }
	} else {
//Form B := alpha*A'*B or B := alpha*conjg(A')*B.
	    if (upper) {
		for (j = 0; j < n; j++) {
		    for (i = m - 1; i >= 0; i--) {
			temp = B[i + j * ldb];
			if (noconj) {
			    if (nounit)
				temp = temp * A[i + i * lda];
			    for (k = 0; k < i; k++) {
				temp = temp + A[k + i * lda] * B[k + j * ldb];
			    }
			} else {
			    if (nounit)
				temp = temp * conj(A[i + i * lda]);
			    for (k = 0; k < i; k++) {
				temp = temp + conj(A[k + i * lda]) * B[k + j * ldb];
			    }
			}
			B[i + j * ldb] = alpha * temp;
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (i = 0; i < m; i++) {
			temp = B[i + j * ldb];
			if (noconj) {
			    if (nounit)
				temp = temp * A[i + i * lda];
			    for (k = i + 1; k < m; k++) {
				temp = temp + A[k + i * lda] * B[k + j * ldb];
			    }
			} else {
			    if (nounit)
				temp = temp * conj(A[i + i * lda]);
			    for (k = i + 1; k < m; k++) {
				temp = temp + conj(A[k + i * lda]) * B[k + j * ldb];
			    }
			}
			B[i + j * ldb] = alpha * temp;
		    }
		}
	    }
	}
    } else {
	if (Mlsame(transa, "N")) {
//Form B := alpha*B*A.
	    if (upper) {
		for (j = n - 1; j >= 0; j--) {
		    temp = alpha;
		    if (nounit)
			temp = temp * A[j + j * lda];
		    for (i = 0; i < m; i++) {
			B[i + j * ldb] = temp * B[i + j * ldb];
		    }
		    for (k = 0; k < j; k++) {
			if (A[k + j * lda] != Zero) {
			    temp = alpha * A[k + j * lda];
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    temp = alpha;
		    if (nounit)
			temp = temp * A[j + j * lda];
		    for (i = 0; i < m; i++) {
			B[i + j * ldb] = temp * B[i + j * ldb];
		    }
		    for (k = j + 1; k < n; k++) {
			if (A[k + j * lda] != Zero) {
			    temp = alpha * A[k + j * lda];
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		}
	    }
	} else {
//form  B := alpha*B*A' or b := alpha*B*conjg(A').
	    if (upper) {
		for (k = 0; k < n; k++) {
		    for (j = 0; j < k; j++) {
			if (A[j + k * lda] != Zero) {
			    if (noconj)
				temp = alpha * A[j + k * lda];
			    else
				temp = alpha * conj(A[j + k * lda]);
			}
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = B[i + j * ldb] + temp * B[i + k * ldb];
			}
		    }
		    temp = alpha;
		    if (nounit) {
			if (noconj)
			    temp = temp * A[k + k * lda];
			else
			    temp = temp * conj(A[k + k * lda]);
		    }
		    if (temp != One) {
			for (i = 0; i < m; i++) {
			    B[i + k * ldb] = temp * B[i + k * ldb];
			}
		    }
		}
	    } else {
		for (k = n - 1; k >= 0; k--) {
		    for (j = k + 1; j < n; j++) {
			if (A[j + k * lda] != Zero) {
			    if (noconj)
				temp = alpha * A[j + k * lda];
			    else
				temp = alpha * conj(A[j + k * lda]);
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] + temp * B[i + k * ldb];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			if (noconj)
			    temp = temp * A[k + k * lda];
			else
			    temp = temp * conj(A[k + k * lda]);
		    }
		    if (temp != One) {
			for (i = 0; i < m; i++) {
			    B[i + k * ldb] = temp * B[i + k * ldb];
			}
		    }
		}
	    }
	}
    }
    return;
}
