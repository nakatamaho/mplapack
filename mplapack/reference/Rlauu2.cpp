/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlauu2.cpp,v 1.9 2010/08/18 07:08:03 nakatamaho Exp $ 
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

void Rlauu2(const char *uplo, INTEGER n, REAL * A, INTEGER lda, INTEGER * info)
{
    REAL aii, One = 1.0;
    INTEGER i, upper;

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
	Mxerbla("Rlauu2", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

    if (upper) {
//Compute the product U * U'.
	for (i = 1; i <= n; i++) {
	    aii = A[i - 1 + (i - 1) * lda];
	    if (i < n) {
		A[(i - 1) + (i - 1) * lda] = Rdot(n - i + 1, &A[(i - 1) + (i - 1) * lda], lda, &A[ (i - 1) + (i - 1) * lda], lda);
		Rgemv("No transpose", i - 1, n - i, One, &A[i * lda], lda, &A[(i - 1) + i * lda], lda, aii, &A[(i - 1) * lda], 1);
	    } else {
		Rscal(i, aii, &A[(i - 1) * lda], 1);
	    }
	}
    } else {
//Compute the product L' * L.
	for (i = 1; i <= n; i++) {
	    aii = A[(i - 1) + (i - 1) * lda];
	    if (i < n) {
		A[(i - 1) + (i - 1) * lda] = Rdot(n - i + 1, &A[(i - 1) + (i - 1) * lda], 1, &A[(i - 1) + (i - 1) * lda], 1);
		Rgemv("Transpose", n - i, i - 1, One, &A[i], lda, &A[i + (i - 1) * lda], 1, aii, &A[i - 1], lda);
	    } else {
		Rscal(i, aii, &A[(i - 1)], lda);
	    }
	}
    }
    return;
}
