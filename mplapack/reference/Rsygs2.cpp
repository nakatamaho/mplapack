/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsygs2.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsygs2(INTEGER itype, const char *uplo, INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, INTEGER * info)
{
    INTEGER k;
    REAL ct, akk, bkk;
    INTEGER upper;
    REAL One = 1.0, Half = 0.5;
    REAL mtemp1;

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
	Mxerbla("Rsygs2", -(*info));
	return;
    }

    if (itype == 1) {
	if (upper) {
//Compute inv(U')*A*inv(U)
	    for (k = 0; k < n; k++) {
//Update the upper triangle of A(k:n,k:n)
		akk = A[k + k * lda];
		bkk = B[k + k * ldb];
		akk /= bkk * bkk;
		A[k + k * lda] = akk;
		if (k < n) {
		    mtemp1 = One / bkk;
		    Rscal(n - k, mtemp1, &A[k + (k + 1) * lda], lda);
		    ct = akk * -Half;
		    Raxpy(n - k, ct, &B[k + (k + 1) * ldb], ldb, &A[k + (k + 1) * lda], lda);
		    Rsyr2(uplo, n - k, -One, &A[k + (k + 1) * lda], lda, &B[k + (k + 1) * ldb], ldb, &A[k + 1 + (k + 1)
													* lda], lda);

		    Raxpy(n - k, ct, &B[k + (k + 1) * ldb], ldb, &A[k + (k + 1) * lda], lda);
		    Rtrsv(uplo, "Transpose", "Non-unit", n - k, &B[k + 1 + (k + 1) * ldb], ldb, &A[k + (k + 1) * lda], lda);
		}

	    }
	} else {
//Compute inv(L)*A*inv(L') 
	    for (k = 0; k < n; k++) {
//Update the lower triangle of A(k:n,k:n)
		akk = A[k + k * lda];
		bkk = B[k + k * ldb];
		akk /= bkk * bkk;
		A[k + k * lda] = akk;
		if (k < n) {
		    mtemp1 = One / bkk;
		    Rscal(n - k, mtemp1, &A[k + 1 + k * lda], 1);
		    ct = akk * -Half;
		    Raxpy(n - k, ct, &B[k + 1 + k * ldb], 1, &A[k + 1 + k * lda], 1);
		    Rsyr2(uplo, n - k, -One, &A[k + 1 + k * lda], 1, &B[k + 1 + k * ldb], 1, &A[k + 1 + (k + 1)
												* lda], lda);
		    Raxpy(n - k, ct, &B[k + 1 + k * ldb], 1, &A[k + 1 + k * lda], 1);
		    Rtrsv(uplo, "No transpose", "Non-unit", n - k, &B[k + 1 + (k + 1) * ldb], ldb, &A[k + 1 + k * lda], 1);
		}
	    }
	}
    } else {
	if (upper) {
//Compute U*A*U'
	    for (k = 0; k < n; k++) {
//Update the upper triangle of A(1:k,1:k)
		akk = A[k + k * lda];
		bkk = B[k + k * ldb];
		Rtrmv(uplo, "No transpose", "Non-unit", k - 1, &B[0], ldb, &A[k * lda], 1);
		ct = akk * Half;
		Raxpy(k - 1, ct, &B[k * ldb + 1], 1, &A[k * lda], 1);
		Rsyr2(uplo, k - 1, One, &A[k * lda], 1, &B[k * ldb + 1], 1, &A[0], lda);
		Raxpy(k - 1, ct, &B[k * ldb + 1], 1, &A[k * lda], 1);
		Rscal(k - 1, bkk, &A[k * lda], 1);
		A[k + k * lda] = akk * (bkk * bkk);

	    }
	} else {
//Compute L'*A*L
	    for (k = 0; k < n; k++) {
//Update the lower triangle of A(1:k,1:k)
		akk = A[k + k * lda];
		bkk = B[k + k * ldb];
		Rtrmv(uplo, "Transpose", "Non-unit", k - 1, &B[0], ldb, &A[k + lda], lda);
		ct = akk * Half;
		Raxpy(k - 1, ct, &B[k + ldb], ldb, &A[k + lda], lda);
		Rsyr2(uplo, k - 1, One, &A[k + lda], lda, &B[k + ldb], ldb, &A[0], lda);
		Raxpy(k - 1, ct, &B[k + ldb], ldb, &A[k + lda], lda);
		Rscal(k - 1, bkk, &A[k + lda], lda);
		A[k + k * lda] = akk * (bkk * bkk);

	    }
	}
    }
    return;
}
