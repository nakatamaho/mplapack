/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Checon.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Checon(const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * ipiv, REAL anorm, REAL * rcond, COMPLEX * work, INTEGER * info)
{
    INTEGER i, kase;
    INTEGER isave[3];
    INTEGER upper;
    REAL ainvnm;
    REAL Zero = 0.0, One = 1.0;
//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    } else if (anorm < Zero) {
	*info = -6;
    }
    if (*info != 0) {
	Mxerbla("Checon", -(*info));
	return;
    }
//Quick return if possible
    *rcond = Zero;
    if (n == 0) {
	*rcond = One;
	return;
    } else if (anorm <= Zero) {
	return;
    }
//Check that the diagonal matrix D is nonsingular.
    if (upper) {
//Upper triangular storage: examine D from bottom to top
	for (i = n; i >= 1; i--) {
	    if (ipiv[i] > 0 && A[i + i * lda] == Zero) {
		return;
	    }
	}
    } else {
//Lower triangular storage: examine D from top to bottom.
	for (i = 0; i < n; i++) {
	    if (ipiv[i] > 0 && A[i + i * lda] == Zero) {
		return;
	    }
	}
    }
//Estimate the 1-norm of the inverse.
    kase = 0;
  L30:
    Clacn2(n, &work[n + 1], &work[0], &ainvnm, &kase, isave);
    if (kase != 0) {
//Multiply by inv(L*D*L') or inv(U*D*U').
	Chetrs(uplo, n, 1, &A[0], lda, &ipiv[1], &work[0], n, info);
	goto L30;
    }
//Compute the estimate of the reciprocal condition number.
    if (ainvnm != Zero) {
	*rcond = One / ainvnm / anorm;
    }
    return;
}
