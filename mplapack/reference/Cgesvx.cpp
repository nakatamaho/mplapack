/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgesvx.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Cgesvx(const char *fact, const char *trans, INTEGER n, INTEGER nrhs,
	    COMPLEX * A, INTEGER lda, COMPLEX * af, INTEGER ldaf, INTEGER * ipiv, char *equed, REAL * r, REAL * c,
	    COMPLEX * B, INTEGER ldb, COMPLEX * x, INTEGER ldx, REAL * rcond, REAL * ferr, REAL * berr, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j;
    REAL amax;
    char norm;
    REAL rcmin, rcmax, anorm;
    LOGICAL equil;
    REAL colcnd;
    LOGICAL nofact;
    REAL bignum = 0.0;
    INTEGER infequ;
    LOGICAL colequ;
    REAL rowcnd;
    LOGICAL notran;
    REAL smlnum = 0.0;
    LOGICAL rowequ;
    REAL rpvgrw;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    nofact = Mlsame(fact, "N");
    equil = Mlsame(fact, "E");
    notran = Mlsame(trans, "N");
    if (nofact || equil) {
	*equed = 'N';
	rowequ = MFALSE;
	colequ = MFALSE;
    } else {
	rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
	colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
	smlnum = Rlamch("Safe minimum");
	bignum = One / smlnum;
    }
//Test the input parameters.
    if (!nofact && !equil && !Mlsame(fact, "F")) {
	*info = -1;
    } else if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (nrhs < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldaf < max((INTEGER) 1, n)) {
	*info = -8;
    } else if (Mlsame(fact, "F") && !(rowequ || colequ || Mlsame(equed, "N"))) {
	*info = -10;
    } else {
	if (rowequ) {
	    rcmin = bignum;
	    rcmax = Zero;
	    for (j = 0; j < n; j++) {
		mtemp1 = rcmin, mtemp2 = r[j];
		rcmin = min(mtemp1, mtemp2);
		mtemp1 = rcmax, mtemp2 = r[j];
		rcmax = max(mtemp1, mtemp2);
	    }
	    if (rcmin <= Zero) {
		*info = -11;
	    } else if (n > 0) {
		rowcnd = max(rcmin, smlnum) / min(rcmax, bignum);
	    } else {
		rowcnd = One;
	    }
	}
	if (colequ && *info == 0) {
	    rcmin = bignum;
	    rcmax = Zero;
	    for (j = 0; j < n; j++) {
		mtemp1 = rcmin, mtemp2 = c[j];
		rcmin = min(mtemp1, mtemp2);
		mtemp1 = rcmax, mtemp2 = c[j];
		rcmax = max(mtemp1, mtemp2);
	    }
	    if (rcmin <= Zero) {
		*info = -12;
	    } else if (n > 0) {
		colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
	    } else {
		colcnd = One;
	    }
	}
	if (*info == 0) {
	    if (ldb < max((INTEGER) 1, n)) {
		*info = -14;
	    } else if (ldx < max((INTEGER) 1, n)) {
		*info = -16;
	    }
	}
    }
    if (*info != 0) {
	Mxerbla("Cgesvx", -(*info));
	return;
    }
    if (equil) {
//Compute row and column scalings to equilibrate the matrix A.
	Cgeequ(n, n, &A[0], lda, &r[1], &c[1], &rowcnd, &colcnd, &amax, &infequ);
	if (infequ == 0) {
//Equilibrate the matrix.
	    Claqge(n, n, &A[0], lda, &r[1], &c[1], rowcnd, colcnd, amax, equed);
	    rowequ = Mlsame(equed, "R") || Mlsame(equed, "B");
	    colequ = Mlsame(equed, "C") || Mlsame(equed, "B");
	}
    }
//Scale the right hand side.
    if (notran) {
	if (rowequ) {
	    for (j = 0; j < nrhs; j++) {
		for (i = 0; i < n; i++) {
		    B[i + j * ldb] = r[i] * B[i + j * ldb];
		}
	    }
	}
    } else if (colequ) {
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < n; i++) {
		B[i + j * ldb] = c[i] * B[i + j * ldb];
	    }
	}
    }
    if (nofact || equil) {
//Compute the LU factorization of A.
	Clacpy("Full", n, n, &A[0], lda, &af[0], ldaf);
	Cgetrf(n, n, &af[0], ldaf, &ipiv[1], info);
//Return if INFO is non-zero.
	if (*info > 0) {
//Compute the reciprocal pivot growth factor of the
//leading rank-deficient INFO columns of A.
	    rpvgrw = Clantr("M", "U", "N", *info, *info, &af[0], ldaf, &rwork[1]);
	    if (rpvgrw == Zero) {
		rpvgrw = One;
	    } else {
		rpvgrw = Clange("M", n, *info, &A[0], lda, &rwork[1]) / rpvgrw;
	    }
	    rwork[1] = rpvgrw;
	    *rcond = Zero;
	    return;
	}
    }
//Compute the norm of the matrix A and the
//reciprocal pivot growth factor RPVGRW.
    if (notran) {
	norm = '1';
    } else {
	norm = 'I';
    }
    anorm = Clange((const char *) norm, n, n, &A[0], lda, &rwork[1]);
    rpvgrw = Clantr("M", "U", "N", n, n, &af[0], ldaf, &rwork[1]);
    if (rpvgrw == Zero) {
	rpvgrw = One;
    } else {
	rpvgrw = Clange("M", n, n, &A[0], lda, &rwork[1]) / rpvgrw;
    }
//Compute the reciprocal of the condition number of A.
    Cgecon((const char *) norm, n, &af[0], ldaf, anorm, rcond, &work[0], &rwork[1], info);
//Compute the solution matrix X.
    Clacpy("Full", n, nrhs, &B[0], ldb, &x[0], ldx);
    Cgetrs(trans, n, nrhs, &af[0], ldaf, &ipiv[1], &x[0], ldx, info);
//Use iterative refinement to improve the computed solution and
//compute error bounds and backward error estimates for it.
    Cgerfs(trans, n, nrhs, &A[0], lda, &af[0], ldaf, &ipiv[1], &B[0], ldb, &x[0], ldx, &ferr[1], &berr[1], &work[0], &rwork[1], info);
//Transform the solution matrix X to a solution of the original
//system.
    if (notran) {
	if (colequ) {
	    for (j = 0; j < nrhs; j++) {
		for (i = 0; i < n; i++) {
		    x[i + j * ldx] = c[i] * x[i + j * ldx];
		}
	    }
	    for (j = 0; j < nrhs; j++) {
		ferr[j] = ferr[j] / colcnd;
	    }
	}
    } else if (rowequ) {
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < n; i++) {
		x[i + j * ldx] = r[i] * x[i + j * ldx];
	    }
	}
	for (j = 0; j < nrhs; j++) {
	    ferr[j] = ferr[j] / rowcnd;
	}
    }
//Set INFO = N+1 if the matrix is singular to working precision.
    if (*rcond < Rlamch("Epsilon")) {
	*info = n + 1;
    }
    rwork[1] = rpvgrw;
    return;
}
