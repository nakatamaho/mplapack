/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chegst.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chegst(INTEGER itype, const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, INTEGER * info)
{
    INTEGER k, kb, nb;
    INTEGER upper;
    REAL One = 1.0;
    REAL Half = 0.5;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (itype < 1 || itype > 3) {
	*info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -7;
    }
    if (*info != 0) {
	Mxerbla("Chegst", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Determine the block size for this environment.
    nb = iMlaenv(1, "Chegst", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
//Use unblocked code
	Chegs2(itype, uplo, n, &A[0], lda, &B[0], ldb, info);
    } else {
//Use blocked code
	if (itype == 1) {
	    if (upper) {
//Compute inv(U')*A*inv(U)
		for (k = 0; k <= n; k += nb) {
		    kb = min(n - k + 1, nb);
//Update the upper triangle of A(k:n,k:n)
		    Chegs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		    if (k + kb <= n) {
			Ctrsm("Left", uplo, "Conjugate transpose", "Non-unit", kb, n - k - kb + 1, One, &B[k + k * ldb], ldb, &A[k + (k + kb) * lda], lda);
			Chemm("Left", uplo, kb, n - k - kb + 1, (COMPLEX) - Half, &A[k + k * lda], lda, &B[k + (k + kb) * ldb], ldb, (COMPLEX) One, &A[k + (k + kb) * lda], lda);
			Cher2k(uplo, "Conjugate transpose", n - k - kb + 1, kb, (COMPLEX) - One, &A[k + (k + kb) * lda], lda,
			       &B[k + (k + kb) * ldb], ldb, One, &A[k + kb + (k + kb) * lda], lda);
			Chemm("Left", uplo, kb, n - k - kb + 1, (COMPLEX) - Half, &A[k + k * lda], lda, &B[k + (k + kb) * ldb], ldb, (COMPLEX) One, &A[k + (k + kb) * lda], lda);
			Ctrsm("Right", uplo, "No transpose", "Non-unit", kb, n - k - kb + 1, (COMPLEX) One, &B[k + kb + (k + kb) * ldb], ldb, &A[k + (k + kb) * lda], lda);
		    }
		}
	    } else {
//Compute inv(L)*A*inv(L')
		for (k = 0; k <= n; k += nb) {
		    kb = min(n - k + 1, nb);
// Update the lower triangle of A(k:n,k:n)
		    Chegs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		    if (k + kb <= n) {
			Ctrsm("Right", uplo, "Conjugate transpose", "Non-un" "it", n - k - kb + 1, kb, One, &B[k + k * ldb], ldb, &A[k + kb + k * lda], lda);
			Chemm("Right", uplo, n - k - kb + 1, kb, (COMPLEX) - Half, &A[k + k * lda], lda, &B[k + kb + k * ldb], ldb, One, &A[k + kb + k * lda], lda);
			Cher2k(uplo, "No transpose", n - k - kb + 1, kb, (COMPLEX) - One, &A[k + kb + k * lda], lda,
			       &B[k + kb + k * ldb], ldb, One, &A[k + kb + (k + kb) * lda], lda);
			Chemm("Right", uplo, n - k - kb + 1, kb, (COMPLEX) - Half, &A[k + k * lda], lda, &B[k + kb + k * ldb], ldb, One, &A[k + kb + k * lda], lda);
			Ctrsm("Left", uplo, "No transpose", "Non-unit", n - k - kb + 1, kb, One, &B[k + kb + (k + kb) * ldb], ldb, &A[k + kb + k * lda], lda);
		    }
		}
	    }
	} else {
	    if (upper) {
//Compute U*A*U'
		for (k = 0; k <= n; k = k + nb) {
		    kb = min(n - k + 1, nb);
//Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
		    Ctrmm("Left", uplo, "No transpose", "Non-unit", k - 1, kb, One, &B[0], ldb, &A[k * lda], lda);
		    Chemm("Right", uplo, k - 1, kb, (COMPLEX) Half, &A[k + k * lda], lda, &B[k * ldb + 1], ldb, One, &A[k * lda], lda);
		    Cher2k(uplo, "No transpose", k - 1, kb, (COMPLEX) One, &A[k * lda], lda, &B[k * ldb + 1], ldb, One, &A[0], lda);
		    Chemm("Right", uplo, k - 1, kb, (COMPLEX) Half, &A[k + k * lda], lda, &B[k * ldb + 1], ldb, One, &A[k * lda], lda);
		    Ctrmm("Right", uplo, "Conjugate transpose", "Non-unit", k - 1, kb, One, &B[k + k * ldb], ldb, &A[k * lda], lda);
		    Chegs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		}
	    } else {
//Compute L'*A*L
		for (k = 0; k <= n; k += nb) {
		    kb = min(n - k + 1, nb);
//Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
		    Ctrmm("Right", uplo, "No transpose", "Non-unit", kb, k - 1, One, &B[0], ldb, &A[k + lda], lda);
		    Chemm("Left", uplo, kb, k - 1, (COMPLEX) Half, &A[k + k * lda]
			  , lda, &B[k + ldb], ldb, One, &A[k + lda], lda);
		    Cher2k(uplo, "Conjugate transpose", k - 1, kb, One, &A[k + lda], lda, &B[k + ldb], ldb, One, &A[0], lda);
		    Chemm("Left", uplo, kb, k - 1, Half, &A[k + k * lda]
			  , lda, &B[k + ldb], ldb, One, &A[k + lda], lda);
		    Ctrmm("Left", uplo, "Conjugate transpose", "Non-unit", kb, k - 1, One, &B[k + k * ldb], ldb, &A[k + lda], lda);
		    Chegs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		}
	    }
	}
    }
    return;
}
