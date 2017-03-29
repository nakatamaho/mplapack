/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsysvx.cpp,v 1.3 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void
Rsysvx(const char *fact, const char *uplo, INTEGER n, INTEGER nrhs, REAL * A, INTEGER lda, REAL * af, INTEGER ldaf, INTEGER * ipiv,
       REAL * B, INTEGER ldb, REAL * x, INTEGER ldx, REAL * rcond, REAL * ferr, REAL * berr, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER nb;
    REAL anorm;
    INTEGER nofact;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL Zero = 0.0;

//Test the input parameters.
    *info = 0;
    nofact = Mlsame(fact, "N");
    lquery = lwork == -1;
    if (!nofact && !Mlsame(fact, "F")) {
	*info = -1;
    } else if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (nrhs < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldaf < max((INTEGER) 1, n)) {
	*info = -8;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -11;
    } else if (ldx < max((INTEGER) 1, n)) {
	*info = -13;
    } else {
	if (lwork < max((INTEGER) 1, n * 3) && !lquery) {
	    *info = -18;
	}
    }
    if (*info == 0) {
	lwkopt = max((INTEGER) 1, n * 3);
	if (nofact) {
	    nb = iMlaenv(1, "Rsytrf", uplo, n, -1, -1, -1);
	    lwkopt = max(lwkopt, n * nb);
	}
	work[1] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Rsysvx", -(*info));
	return;
    } else if (lquery) {
	return;
    }
    if (nofact) {
//Compute the factorization A = U*D*U' or A = L*D*L'.
	Rlacpy(uplo, n, n, &A[0], lda, &af[0], ldaf);
	Rsytrf(uplo, n, &af[0], ldaf, &ipiv[1], &work[0], lwork, info);
//Return if INFO is non-zero.
	if (*info > 0) {
	    *rcond = Zero;
	    return;
	}
    }
//Compute the norm of the matrix A.
    anorm = Rlansy("I", uplo, n, &A[0], lda, &work[0]);
//Compute the reciprocal of the condition number of A.
    Rsycon(uplo, n, &af[0], ldaf, &ipiv[1], anorm, rcond, &work[0], &iwork[1], info);
//Compute the solution vectors X.
    Rlacpy("Full", n, nrhs, &B[0], ldb, &x[0], ldx);
    Rsytrs(uplo, n, nrhs, &af[0], ldaf, &ipiv[1], &x[0], ldx, info);
//Use iterative refinement to improve the computed solutions and
//compute error bounds and backward error estimates for them.
    Rsyrfs(uplo, n, nrhs, &A[0], lda, &af[0], ldaf, &ipiv[1], &B[0], ldb, &x[0], ldx, &ferr[1], &berr[1], &work[0]
	   , &iwork[1], info);
//Set INFO = N+1 if the matrix is singular to working precision.
    if (*rcond < Rlamch("Epsilon")) {
	*info = n + 1;
    }
    work[1] = lwkopt;
    return;
}
