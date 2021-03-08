/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsygst.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsygst(INTEGER itype, const char *uplo, INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, INTEGER * info)
{
    INTEGER k, kb, nb;
    REAL One = 1.0, Half = 0.5;
    INTEGER upper;

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
	Mxerbla("Rsygst", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

//Determine the block size for this environment.
    nb = iMlaenv(1, "Rsygst", uplo, n, -1, -1, -1);

    if (nb <= 1 || nb >= n) {
//Use unblocked code
	Rsygs2(itype, uplo, n, &A[0], lda, &B[0], ldb, info);

    } else {
//Use blocked code
	if (itype == 1) {
	    if (upper) {
//Compute inv(U')*A*inv(U)
		for (k = 0; k < n; k += nb) {
		    kb = min(n - k + 1, nb);
//Update the upper triangle of A(k:n,k:n)

		    Rsygs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		    if (k + kb <= n) {

			Rtrsm("Left", uplo, "Transpose", "Non-unit", kb, n - k - kb + 1, One, &B[k + k * ldb], ldb, &A[k + (k + kb) * lda], lda);

			Rsymm("Left", uplo, kb, n - k - kb + 1, -Half, &A[k + k * lda], lda, &B[k + (k + kb) * ldb], ldb, One, &A[k + (k + kb) * lda], lda);
			Rsyr2k(uplo, "Transpose", n - k - kb + 1, kb, -One, &A[k + (k + kb) * lda], lda, &B[k + (k + kb) * ldb], ldb, One, &A[k + kb + (k + kb) * lda], lda);

			Rsymm("Left", uplo, kb, n - k - kb + 1, -Half, &A[k + k * lda], lda, &B[k + (k + kb) * ldb], ldb, One, &A[k + (k + kb) * lda], lda);

			Rtrsm("Right", uplo, "No transpose", "Non-unit", kb, n - k - kb + 1, One, &B[k + kb + (k + kb) * ldb]
			      , ldb, &A[k + (k + kb) * lda], lda);
		    }

		}
	    } else {
//Compute inv(L)*A*inv(L')
		for (k = 0; k < n; k += nb) {
		    kb = min(n - k + 1, nb);
//Update the lower triangle of A(k:n,k:n)
		    Rsygs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		    if (k + kb <= n) {
			Rtrsm("Right", uplo, "Transpose", "Non-unit", n - k - kb + 1, kb, One, &B[k + k * ldb], ldb, &A[k + kb + k * lda], lda);
			Rsymm("Right", uplo, n - k - kb + 1, kb, -Half, &A[k + k * lda], lda, &B[k + kb + k * ldb], ldb, One, &A[k + kb + k * lda], lda);
			Rsyr2k(uplo, "No transpose", n - k - kb + 1, kb, -One, &A[k + kb + k * lda], lda, &B[k + kb + k * ldb], ldb, One, &A[k + kb + (k + kb) * lda], lda);
			Rsymm("Right", uplo, n - k - kb + 1, kb, -Half, &A[k + k * lda], lda, &B[k + kb + k * ldb], ldb, One, &A[k + kb + k * lda], lda);
			Rtrsm("Left", uplo, "No transpose", "Non-unit", n - k - kb + 1, kb, One, &B[k + kb + (k + kb) * ldb], ldb, &A[k + kb + k * lda], lda);
		    }
		}
	    }
	} else {
	    if (upper) {
//Compute U*A*U'
		for (k = 0; k < n; k += nb) {
		    kb = min(n - k + 1, nb);
//Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
		    Rtrmm("Left", uplo, "No transpose", "Non-unit", k - 1, kb, One, &B[0], ldb, &A[k * lda], lda);
		    Rsymm("Right", uplo, k - 1, kb, Half, &A[k + k * lda], lda, &B[k * ldb + 1], ldb, One, &A[k * lda], lda);
		    Rsyr2k(uplo, "No transpose", k - 1, kb, One, &A[k * lda], lda, &B[k * ldb + 1], ldb, One, &A[0], lda);
		    Rsymm("Right", uplo, k - 1, kb, Half, &A[k + k * lda], lda, &B[k * ldb + 1], ldb, One, &A[k * lda], lda);
		    Rtrmm("Right", uplo, "Transpose", "Non-unit", k - 1, kb, One, &B[k + k * ldb], ldb, &A[k * lda], lda);
		    Rsygs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		}
	    } else {
//Compute L'*A*L
		for (k = 0; k < n; k += nb) {
		    kb = min(n - k + 1, nb);
//Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
		    Rtrmm("Right", uplo, "No transpose", "Non-unit", kb, k - 1, One, &B[0], ldb, &A[k + lda], lda);
		    Rsymm("Left", uplo, kb, k - 1, Half, &A[k + k * lda], lda, &B[k + ldb], ldb, One, &A[k + lda], lda);
		    Rsyr2k(uplo, "Transpose", k - 1, kb, One, &A[k + lda], lda, &B[k + ldb], ldb, One, &A[0], lda);
		    Rsymm("Left", uplo, kb, k - 1, Half, &A[k + k * lda], lda, &B[k + ldb], ldb, One, &A[k + lda], lda);
		    Rtrmm("Left", uplo, "Transpose", "Non-unit", kb, k - 1, One, &B[k + k * ldb], ldb, &A[k + lda], lda);
		    Rsygs2(itype, uplo, kb, &A[k + k * lda], lda, &B[k + k * ldb], ldb, info);
		}
	    }
	}
    }
    return;
}
