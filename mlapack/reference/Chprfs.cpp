/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chprfs.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chprfs(const char *uplo, INTEGER n, INTEGER nrhs,
	    COMPLEX * ap, COMPLEX * afp, INTEGER * ipiv, COMPLEX * B, INTEGER ldb, COMPLEX * x, INTEGER ldx, REAL * ferr, REAL * berr, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j, k;
    REAL s;
    INTEGER ik, kk;
    REAL xk;
    INTEGER nz;
    REAL eps;
    INTEGER kase;
    REAL safe1, safe2;
    INTEGER isave[3], count;
    INTEGER upper;
    REAL safmin;
    REAL lstres;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (nrhs < 0) {
	*info = -3;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    } else if (ldx < max((INTEGER) 1, n)) {
	*info = -10;
    }
    if (*info != 0) {
	Mxerbla("Chprfs", -(*info));
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
    nz = n + 1;
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
	Ccopy(n, &B[j * ldb + 1], 1, &work[0], 1);
	Chpmv(uplo, n, (COMPLEX) - One, &ap[1], &x[j * ldx + 1], 1, (COMPLEX) One, &work[0], 1);
//Compute componentwise relative backward error from formula
//max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
//where abs(Z) is the componentwise absolute value of the matrix
//or vector Z.  If the i-th component of the denominator is less
//than SAFE2, then SAFE1 is added to the i-th components of the
//numerator and denominator before dividing.
	for (i = 0; i < n; i++) {
	    rwork[i] = abs(B[i + j * ldb].real()) + abs(B[i + j * ldb].imag());
	}
//Compute abs(A)*abs(X) + abs(B).
	kk = 0;
	if (upper) {
	    for (k = 0; k < n; k++) {
		s = Zero;
		xk = abs(x[k + j * ldx].real()) + abs(x[k + j * ldx].imag());
		ik = kk;
		for (i = 0; i < k - 1; i++) {
		    rwork[i] = rwork[i] + abs(ap[ik].real()) + abs(ap[ik].imag());
		    s = s + abs(ap[ik].real()) + abs(ap[ik].imag()) * (abs(x[i + j * ldx].real()) + abs(x[i + j * ldx].imag()));
		    ik++;
		}
		rwork[k] = rwork[k] + abs(ap[kk + k - 1].real()) * xk + s;
		kk = kk + k;
	    }
	} else {
	    for (k = 0; k < n; k++) {
		s = Zero;
		xk = abs(x[k + j * ldx].real()) + abs(x[k + j * ldx].imag());
		rwork[k] = rwork[k] + abs(ap[kk].real()) * xk;
		ik = kk + 1;
		for (i = k + 1; i <= n; i++) {
		    rwork[i] = rwork[i] + (abs(ap[ik].real()) + abs(ap[ik].imag())) * xk;
		    s = s + (abs(ap[ik].real()) + abs(ap[ik].imag())) * (abs(x[i + j * ldx].real()) + abs(x[i + j * ldx].imag()));
		    ik++;
		}
		rwork[k] = rwork[k] + s;
		kk = kk + n - k + 1;
	    }
	}
	s = Zero;
	for (i = 0; i < n; i++) {
	    if (rwork[i] > safe2) {
		mtemp1 = s, mtemp2 = (abs(work[i].real()) + abs(work[i].imag())) / rwork[i];
		s = max(mtemp1, mtemp2);
	    } else {
		mtemp1 = s, mtemp2 = (abs(work[i].real()) + abs(work[i].imag()) + safe1) / (rwork[i] + safe1);
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
	    Chptrs(uplo, n, 1, &afp[1], &ipiv[1], &work[0], n, info);
	    Caxpy(n, One, &work[0], 1, &x[j * ldx + 1], 1);
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
//Bound error from formula
//norm(X - XTRUE) / norm(X) .le. FERR =
//norm( abs(inv(A))*
//   ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)
//where
//  norm(Z) is the magnitude of the largest component of Z
//  inv(A) is the inverse of A
//  abs(Z) is the componentwise absolute value of the matrix or
//     vector Z
//  NZ is the maximum number of nonzeros in any row of A, plus 1
//  EPS is machine epsilon
//The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
//is incremented by SAFE1 if the i-th component of
//abs(A)*abs(X) + abs(B) is less than SAFETwo
//Use ZLACN2 to estimate the infinity-norm of the matrix
//   inv(A) * diag(W),
//where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
	for (i = 0; i < n; i++) {
	    if (rwork[i] > safe2) {
		rwork[i] = abs(work[i].real()) + abs(work[i].imag()) + nz * eps * rwork[i];
	    } else {
		rwork[i] = abs(work[i].real()) + abs(work[i].imag()) + nz * eps * rwork[i] + safe1;
	    }
	}
	kase = 0;
      L100:
	Clacn2(n, &work[n + 1], &work[0], &ferr[j], &kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
//Multiply by diag(W)*inv(A').
		Chptrs(uplo, n, 1, &afp[1], &ipiv[1], &work[0], n, info);
		for (i = 0; i < n; i++) {
		    work[i] = rwork[i] * work[i];
		}
	    } else if (kase == 2) {
//Multiply by inv(A)*diag(W).
		for (i = 0; i < n; i++) {
		    work[i] = rwork[i] * work[i];
		}
		Chptrs(uplo, n, 1, &afp[1], &ipiv[1], &work[0], n, info);
	    }
	    goto L100;
	}
//Normalize error.
	lstres = Zero;
	for (i = 0; i < n; i++) {
	    mtemp1 = lstres, mtemp2 = abs(x[i + j * ldx].real()) + abs(x[i + j * ldx].imag());
	    lstres = max(mtemp1, mtemp2);
	}
	if (lstres != Zero) {
	    ferr[j] = ferr[j] / lstres;
	}
    }
    return;
}
