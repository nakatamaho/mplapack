/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rposvx.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rposvx(const char *fact, const char *uplo, INTEGER n, INTEGER nrhs, REAL * A, INTEGER lda, REAL * af, INTEGER ldaf,
	    char *equed, REAL * s, REAL * B, INTEGER ldb, REAL * x, INTEGER ldx, REAL * rcond, REAL * ferr, REAL * berr, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j;
    REAL amax, smin, smax;
    REAL scond, anorm;
    LOGICAL equil, rcequ;
    LOGICAL nofact;
    REAL bignum = 0.0;
    INTEGER infequ;
    REAL smlnum = 0.0;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    nofact = Mlsame(fact, "N");
    equil = Mlsame(fact, "E");
    if (nofact || equil) {
	*equed = 'N';
	rcequ = MFALSE;
    } else {
	rcequ = Mlsame(equed, "Y");
	smlnum = Rlamch("Safe minimum");
	bignum = One / smlnum;
    }
//Test the input parameters.
    if (!nofact && !equil && !Mlsame(fact, "F")) {
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
    } else if (Mlsame(fact, "F") && !(rcequ || Mlsame(equed, "N"))) {
	*info = -9;
    } else {
	if (rcequ) {
	    smin = bignum;
	    smax = Zero;
	    for (j = 0; j < n; j++) {
		mtemp1 = smin, mtemp2 = s[j];
		smin = min(mtemp1, mtemp2);
		mtemp1 = smax, mtemp2 = s[j];
		smax = max(mtemp1, mtemp2);
	    }
	    if (smin <= Zero) {
		*info = -10;
	    } else if (n > 0) {
		scond = max(smin, smlnum) / min(smax, bignum);
	    } else {
		scond = One;
	    }
	}
	if (*info == 0) {
	    if (ldb < max((INTEGER) 1, n)) {
		*info = -12;
	    } else if (ldx < max((INTEGER) 1, n)) {
		*info = -14;
	    }
	}
    }
    if (*info != 0) {
	Mxerbla("Rposvx", -(*info));
	return;
    }
    if (equil) {
//Compute row and column scalings to equilibrate the matrix A.
	Rpoequ(n, &A[0], lda, &s[1], &scond, &amax, &infequ);
	if (infequ == 0) {
//Equilibrate the matrix.
	    Rlaqsy(uplo, n, &A[0], lda, &s[1], scond, amax, equed);
	    rcequ = Mlsame(equed, "Y");
	}
    }
//Scale the right hand side.
    if (rcequ) {
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < n; i++) {
		B[i + j * ldb] = s[i] * B[i + j * ldb];
	    }
	}
    }
    if (nofact || equil) {
//Compute the Cholesky factorization A = U'*U or A = L*L'.
	Rlacpy(uplo, n, n, &A[0], lda, &af[0], ldaf);
	Rpotrf(uplo, n, &af[0], ldaf, info);
//Return if INFO is non-zero.
	if (*info > 0) {
	    *rcond = Zero;
	    return;
	}
    }
//Compute the norm of the matrix A.
    anorm = Rlansy("1", uplo, n, &A[0], lda, &work[0]);
//Compute the reciprocal of the condition number of A.
    Rpocon(uplo, n, &af[0], ldaf, anorm, rcond, &work[0], &iwork[1], info);
//Compute the solution matrix X.
    Rlacpy("Full", n, nrhs, &B[0], ldb, &x[0], ldx);
    Rpotrs(uplo, n, nrhs, &af[0], ldaf, &x[0], ldx, info);
//Use iterative refinement to improve the computed solution and
//compute error bounds and backward error estimates for it.
    Rporfs(uplo, n, nrhs, &A[0], lda, &af[0], ldaf, &B[0], ldb, &x[0], ldx, &ferr[1], &berr[1], &work[0], &iwork[1], info);
//Transform the solution matrix X to a solution of the original
//system.
    if (rcequ) {
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < n; i++) {
		x[i + j * ldx] = s[i] * x[i + j * ldx];

	    }
	}
	for (j = 0; j < nrhs; j++) {
	    ferr[j] = ferr[j] / scond;
	}
    }
//Set INFO = N+1 if the matrix is singular to working precision.
    if (*rcond < Rlamch("Epsilon")) {
	*info = n + 1;
    }
    return;
}
