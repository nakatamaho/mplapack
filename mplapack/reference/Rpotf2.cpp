/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rpotf2.cpp,v 1.10 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rpotf2(const char *uplo, INTEGER n, REAL * A, INTEGER lda, INTEGER * info)
{
    INTEGER j, upper, success = 1;
    REAL ajj;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rpotf2", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (upper) {
//Compute the Cholesky factorization A = U'*U.
	for (j = 0; j < n; j++) {
//Compute U(J,J) and test for non-positive-definiteness.
	    ajj = A[j + j * lda] - Rdot(j, &A[j * lda], 1, &A[j * lda], 1);
	    if (ajj <= Zero) {
		A[j + j * lda] = ajj;
		success = 0;
		break;
	    }
	    ajj = sqrt(ajj);
	    A[j + j * lda] = ajj;
//Compute elements J+1:N of row J.
	    if (j < n) {
		Rgemv("Transpose", j, n - j - 1, -One, &A[(j + 1) * lda], lda, &A[j * lda], 1, One, &A[j + (j + 1) * lda], lda);
		Rscal(n - j - 1, One / ajj, &A[j + (j + 1) * lda], lda);
	    }
	}
    } else {
//Compute the Cholesky factorization A = L*L'.
	for (j = 0; j < n; j++) {
//Compute L(J,J) and test for non-positive-definiteness.
	    ajj = A[j + j * lda] - Rdot(j, &A[j], lda, &A[j], lda);
	    if (ajj <= Zero) {
		A[j + j * lda] = ajj;
		success = 0;
		break;
	    }
	    ajj = sqrt(ajj);
	    A[j + j * lda] = ajj;
//Compute elements J+1:N of column J.
	    if (j < n) {
		Rgemv("No transpose", n - j - 1, j, -One, &A[j + 1], lda, &A[j], lda, One, &A[j + 1 + j * lda], 1);
		Rscal(n - j - 1, One / ajj, &A[j + 1 + j * lda], 1);
	    }
	}
    }
    if (!success)
	*info = j + 1;
    return;
}
