/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rpbrfs.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rpbrfs(const char *uplo, INTEGER n, INTEGER kd, INTEGER nrhs, REAL * AB, INTEGER ldab, REAL * afb, INTEGER ldafb,
	    REAL * B, INTEGER ldb, REAL * x, INTEGER ldx, REAL * ferr, REAL * berr, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, k, l;
    REAL s, xk;
    INTEGER nz;
    REAL eps;
    INTEGER kase;
    REAL safe1, safe2;
    INTEGER isave[3];
    INTEGER count;
    INTEGER upper;
    REAL safmin;
    REAL lstres;
    REAL mtemp1, mtemp2;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kd < 0) {
	*info = -3;
    } else if (nrhs < 0) {
	*info = -4;
    } else if (ldab < kd + 1) {
	*info = -6;
    } else if (ldafb < kd + 1) {
	*info = -8;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -10;
    } else if (ldx < max((INTEGER) 1, n)) {
	*info = -12;
    }
    if (*info != 0) {
	Mxerbla("Rpbrfs", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0 || nrhs == 0) {
	for (j = 0; j < nrhs; j++) {
	    ferr[j] = Zero;
	    berr[j] = Zero;
	}
	return;
    }
//NZ = maximum number of nonzero elements in each row of A, plus 1
/* Computing MIN */
    nz = min(n + 1, (kd * 2) + 2);
    eps = Rlamch("Epsilon");
    safmin = Rlamch("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
//Do for each right hand side
    for (j = 0; j < nrhs; j++) {
	count = 1;
	lstres = 3.;
      L20:
//Loop until stopping criterion is satisfied.
//Compute residual R = B - A * X
	Rcopy(n, &B[j * ldb + 1], 1, &work[n + 1], 1);
	Rsbmv(uplo, n, kd, -One, &AB[0], ldab, &x[j * ldx + 1], 1, One, &work[n + 1], 1);
//Compute componentwise relative backward error from formula
//max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

//where abs(Z) is the componentwise absolute value of the matrix
//or vector Z.  If the i-th component of the denominator is less
//than SAFE2, then SAFE1 is added to the i-th components of the
//numerator and denominator before dividing.
	for (i = 0; i < n; i++) {
	    work[i] = abs(B[i + j * ldb]);
	}
//Compute abs(A)*abs(X) + abs(B).
	if (upper) {
	    for (k = 0; k < n; k++) {
		s = Zero;
		xk = abs(x[k + j * ldx]);
		l = kd + 1 - k;
		for (i = max((INTEGER) 1, k - kd); i <= k - 1; i++) {
		    work[i] = work[i] + abs(AB[l + i + k * ldab]) * xk;
		    s = s + abs(AB[l + i + k * ldab]) * abs(x[i + j * ldx]);
		}
		work[k] = work[k] + abs(AB[kd + 1 + k * ldab]) * xk + s;
	    }
	} else {
	    for (k = 0; k < n; k++) {
		s = Zero;
		xk = abs(x[k + j * ldx]);
		work[k] = work[k] + abs(AB[k * ldab + 1]) * xk;
		l = 0 - k;
		for (i = k + 1; i <= min(n, k + kd); i++) {
		    work[i] = work[i] + abs(AB[l + i + k * ldab]) * xk;
		    s = s + abs(AB[l + i + k * ldab]) * abs(x[i + j * ldx]);
		}
		work[k] = work[k] + s;
	    }
	}
	s = Zero;
	for (i = 0; i < n; i++) {
	    if (work[i] > safe2) {
		mtemp1 = s, mtemp2 = abs(work[n + i]) / work[i];
		s = max(mtemp1, mtemp2);
	    } else {
		mtemp1 = s, mtemp2 = (abs(work[n + i]) + safe1) / (work[i] + safe1);
		s = max(mtemp1, mtemp2);
	    }
	}
	berr[j] = s;
//Test stopping criterion. Continue iterating if
//   1) The residual BERR(J) is larger than machine epsilon, and
//   2) BERR(J) decreased by at least a factor of 2 during the
//      last iteration, and
//   3) At most ITMAX iterations tried.
	if (berr[j] > eps && berr[j] * Two <= lstres && count <= 5) {
//Update solution and try again.
	    Rpbtrs(uplo, n, kd, 1, &afb[0], ldafb, &work[n + 1], n, info);
	    Raxpy(n, One, &work[n + 1], 1, &x[j * ldx + 1], 1);
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
//Bound error from formula */
//norm(X - XTRUE) / norm(X) .le. FERR = */
//norm( abs(inv(A))* */
//   ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X) */
//where */
//  norm(Z) is the magnitude of the largest component of Z */
//  inv(A) is the inverse of A */
//  abs(Z) is the componentwise absolute value of the matrix or */
//     vector Z */
//  NZ is the maximum number of nonzeros in any row of A, plus 1 */
//  EPS is machine epsilon */
//The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B)) */
//is incremented by SAFE1 if the i-th component of */
//abs(A)*abs(X) + abs(B) is less than SAFE2
//Use DLACN2 to estimate the infinity-norm of the matrix */
//   inv(A) * diag(W), */
//where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */
	for (i = 0; i < n; i++) {
	    if (work[i] > safe2) {
		work[i] = abs(work[n + i]) + nz * eps * work[i];
	    } else {
		work[i] = abs(work[n + i]) + nz * eps * work[i] + safe1;
	    }
	}
	kase = 0;
      L100:
	Rlacn2(n, &work[(n * 2) + 1], &work[n + 1], &iwork[1], &ferr[j], &kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
//Multiply by diag(W)*inv(A').
		Rpbtrs(uplo, n, kd, 1, &afb[0], ldafb, &work[n + 1], n, info);
		for (i = 0; i < n; i++) {
		    work[n + i] = work[n + i] * work[i];
		}
	    } else if (kase == 2) {
//Multiply by inv(A)*diag(W).
		for (i = 0; i < n; i++) {
		    work[n + i] = work[n + i] * work[i];
		}
		Rpbtrs(uplo, n, kd, 1, &afb[0], ldafb, &work[n + 1], n, info);
	    }
	    goto L100;
	}
//Normalize error.
	lstres = Zero;
	for (i = 0; i < n; i++) {
	    mtemp1 = lstres, mtemp2 = abs(x[i + j * ldx]);
	    lstres = max(mtemp1, mtemp2);
	}
	if (lstres != Zero) {
	    ferr[j] = ferr[j] / lstres;
	}
    }
    return;
}
