/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctrtri.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctrtri(const char *uplo, const char *diag, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * info)
{
    INTEGER j, jb, nb, nn;
    INTEGER upper;
    INTEGER nounit;
    char uplo_diag[3];
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    nounit = Mlsame(diag, "N");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (!nounit && !Mlsame(diag, "U")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    }
    if (*info != 0) {
	Mxerbla("Ctrtri", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Check for singularity if non-unit.
    if (nounit) {
	for (*info = 0; *info < n; ++(*info)) {
	    if (A[*info + *info * lda] == Zero) {
		return;
	    }
	}
	*info = 0;
    }
//Determine the block size for this environment.
    uplo_diag[0] = uplo[0];
    uplo_diag[1] = diag[0];
    uplo_diag[2] = '\0';
    nb = iMlaenv(1, "Ctrtri", uplo_diag, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
//Use unblocked code
	Ctrti2(uplo, diag, n, A, lda, info);
    } else {
//Use blocked code
	if (upper) {
//Compute inverse of upper triangular matrix
	    for (j = 1; j <= n; j = j + nb) {
		jb = min(nb, n - j + 1);
//Compute rows 1:j-1 of current block column
		Ctrmm("Left", "Upper", "No transpose", diag, j - 1, jb, (COMPLEX) One, A, lda, &A[0 + (j - 1) * lda], lda);
		Ctrsm("Right", "Upper", "No transpose", diag, j - 1, jb, (COMPLEX) - One, &A[(j - 1) + (j - 1) * lda], lda, &A[0 + (j - 1) * lda], lda);
//Compute inverse of current diagonal block
		Ctrti2("Upper", diag, jb, &A[(j - 1) + (j - 1) * lda], lda, info);
	    }
	} else {
//Compute inverse of lower triangular matrix
	    nn = (n - 1) / nb * nb + 1;
	    for (j = nn; j >= 1; j = j - nb) {
		jb = min(nb, n - j + 1);
		if (j + jb <= n) {
//Compute rows j+jb:n of current block column
		    Ctrmm("Left", "Lower", "No transpose", diag,
			  n - j - jb + 1, jb, (COMPLEX) One, &A[(j + jb - 1) + (j + jb - 1) * lda], lda, &A[(j + jb - 1) + (j - 1) * lda], lda);
		    Ctrsm("Right", "Lower", "No transpose", diag, n - j - jb + 1, jb, (COMPLEX) - One, &A[(j - 1) + (j - 1) * lda], lda, &A[(j + jb - 1) + (j - 1) * lda], lda);
		}
//Compute inverse of current diagonal block
		Ctrti2("Lower", diag, jb, &A[(j - 1) + (j - 1) * lda], lda, info);
	    }
	}
    }
    return;
}
