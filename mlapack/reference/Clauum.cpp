/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clauum.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clauum(const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * info)
{
    INTEGER i, ib, nb;
    REAL One = 1.0;

    *info = 0;
    INTEGER upper = Mlsame(uplo, "U");

    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Clauum", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Determine the block size for this environment.
    nb = iMlaenv(1, "Clauum", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
//Use unblocked code
	Clauu2(uplo, n, &A[0], lda, info);
    } else {
//Use blocked code
	if (upper) {
//Compute the product U * U'.
	    for (i = 1; i <= n; i = i + nb) {
		ib = min(nb, n - i + 1);
		Ctrmm("Right", "Upper", "Conjugate transpose", "Non-unit", i - 1, ib, One, &A[i + i * lda], lda, &A[i * lda], lda);
		Clauu2("Upper", ib, &A[i + i * lda], lda, info);
		if (i + ib <= n) {
		    Cgemm("No transpose", "Conjugate transpose", i - 1, ib, n - i - ib + 1, One, &A[(i + ib) * lda], lda, &A[i + (i + ib) * lda], lda, One, &A[i * lda], lda);
		    Cherk("Upper", "No transpose", ib, n - i - ib + 1, One, &A[i + (i + ib) * lda], lda, One, &A[i + i * lda], lda);
		}
	    }
	} else {
//Compute the product L' * L.
	    for (i = 1; i <= n; i = i + nb) {
		ib = min(nb, n - i + 1);
		Ctrmm("Left", "Lower", "Conjugate transpose", "Non-unit", ib, i - 1, One, &A[i + i * lda], lda, &A[i + lda], lda);
		Clauu2("Lower", ib, &A[i + i * lda], lda, info);
		if (i + ib <= n) {
		    Cgemm("Conjugate transpose", "No transpose", ib, i - 1, n - i - ib + 1, One, &A[i + ib + i * lda], lda, &A[i + ib + lda], lda, One, &A[i + lda], lda);
		    Cherk("Lower", "Conjugate transpose", ib, n - i - ib + 1, One, &A[i + ib + i * lda], lda, One, &A[i + i * lda], lda);
		}
	    }
	}
    }
    return;
}
