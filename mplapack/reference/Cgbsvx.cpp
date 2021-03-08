/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgbsvx.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgbsvx(const char *fact, const char *trans, INTEGER n, INTEGER kl,
	    INTEGER ku, INTEGER nrhs, COMPLEX * AB, INTEGER ldab,
	    COMPLEX * afb, INTEGER ldafb, INTEGER * ipiv, char *equed,
	    REAL * r, REAL * c, COMPLEX * B, INTEGER ldb, COMPLEX * x, INTEGER ldx, REAL * rcond, REAL * ferr, REAL * berr, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j, j1, j2;
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
    COMPLEX mtemp3;
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
    } else if (kl < 0) {
	*info = -4;
    } else if (ku < 0) {
	*info = -5;
    } else if (nrhs < 0) {
	*info = -6;
    } else if (ldab < kl + ku + 1) {
	*info = -8;
    } else if (ldafb < (kl << 1) + ku + 1) {
	*info = -10;
    } else if (Mlsame(fact, "F") && !(rowequ || colequ || Mlsame(equed, "N"))) {
	*info = -12;
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
		*info = -13;
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
		*info = -14;
	    } else if (n > 0) {
		colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
	    } else {
		colcnd = One;
	    }
	}
	if (*info == 0) {
	    if (ldb < max((INTEGER) 1, n)) {
		*info = -16;
	    } else if (ldx < max((INTEGER) 1, n)) {
		*info = -18;
	    }
	}
    }
    if (*info != 0) {
	Mxerbla("Cgbsvx", -(*info));
	return;
    }
    if (equil) {
//Compute row and column scalings to equilibrate the matrix A.
	Cgbequ(n, n, kl, ku, &AB[0], ldab, &r[1], &c[1], &rowcnd, &colcnd, &amax, &infequ);
	if (infequ == 0) {
//Equilibrate the matrix.
	    Claqgb(n, n, kl, ku, &AB[0], ldab, &r[1], &c[1], rowcnd, colcnd, amax, equed);
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
//Compute the LU factorization of the band matrix A.
	for (j = 0; j < n; j++) {
	    j1 = max(j - ku, (INTEGER) 1);
	    j2 = min(j + kl, n);
	    Ccopy(j2 - j1 + 1, &AB[ku + 1 - j + j1 + j * ldab], 1, &afb[kl + ku + 1 - j + j1 + j * ldafb], 1);
	}
	Cgbtrf(n, n, kl, ku, &afb[0], ldafb, &ipiv[1], info);
//Return if INFO is non-zero.
	if (*info > 0) {
//Compute the reciprocal pivot growth factor of the
//leading rank-deficient INFO columns of A.
	    anorm = Zero;
	    for (j = 0; j < *info; j++) {
		for (i = max(ku + 2 - j, (INTEGER) 1); i <= min(n + ku + 1 - j, kl + ku + 1); i++) {
		    mtemp1 = anorm, mtemp2 = abs(AB[i + j * ldab]);
		    anorm = max(mtemp1, mtemp2);
		}
	    }
	    mtemp3 = Clantb("M", "U", "N", *info, min(*info - 1, kl + ku), &afb[max((INTEGER) 1, kl + ku + 2 - *info) + ldafb], ldafb, &rwork[1]);
	    rpvgrw = mtemp3.real();
	    if (rpvgrw == Zero) {
		rpvgrw = One;
	    } else {
		rpvgrw = anorm / rpvgrw;
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
    anorm = Clangb((const char *) norm, n, kl, ku, &AB[0], ldab, &rwork[1]);
    mtemp3 = Clantb("M", "U", "N", n, kl + ku, &afb[0], ldafb, &rwork[1]);
    rpvgrw = mtemp3.real();

    if (rpvgrw == Zero) {
	rpvgrw = One;
    } else {
	rpvgrw = Clangb("M", n, kl, ku, &AB[0], ldab, &rwork[1]) / rpvgrw;
    }
//Compute the reciprocal of the condition number of A.
    Cgbcon((const char *) norm, n, kl, ku, &afb[0], ldafb, &ipiv[1], anorm, rcond, &work[0], &rwork[1], info);
//Compute the solution matrix X.
    Clacpy("Full", n, nrhs, &B[0], ldb, &x[0], ldx);
    Cgbtrs(trans, n, kl, ku, nrhs, &afb[0], ldafb, &ipiv[1], &x[0], ldx, info);
//Use iterative refinement to improve the computed solution and
//compute error bounds and backward error estimates for it.
    Cgbrfs(trans, n, kl, ku, nrhs, &AB[0], ldab, &afb[0], ldafb, &ipiv[1], &B[0], ldb, &x[0], ldx, &ferr[1], &berr[1], &work[0], &rwork[1], info);
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
