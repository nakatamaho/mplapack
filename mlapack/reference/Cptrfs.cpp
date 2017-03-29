/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cptrfs.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cptrfs(const char *uplo, INTEGER n, INTEGER nrhs,
	    REAL * d, COMPLEX * e, REAL * df, COMPLEX * ef,
	    COMPLEX * B, INTEGER ldb, COMPLEX * x, INTEGER ldx, REAL * ferr, REAL * berr, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j;
    REAL s;
    COMPLEX bi, cx, dx, ex;
    INTEGER ix, nz;
    REAL eps, safe1, safe2;
    INTEGER count;
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
	*info = -9;
    } else if (ldx < max((INTEGER) 1, n)) {
	*info = -11;
    }
    if (*info != 0) {
	Mxerbla("Cptrfs", -(*info));
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
    nz = 4;
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
//Compute residual R = B - A * X.  Also compute
//abs(A)*abs(x) + abs(b) for use in the backward error bound.
	if (upper) {
	    if (n == 1) {
		bi = B[j * ldb + 1];
		dx = d[1] * x[j * ldx + 1];
		work[1] = bi - dx;
		rwork[1] = Cabs1(bi) + Cabs1(dx);
	    } else {
		bi = B[j * ldb + 1];
		dx = d[1] * x[j * ldx + 1];
		ex = e[1] * x[j * ldx + 2];
		work[1] = bi - dx - ex;
		rwork[1] = Cabs1(bi) + Cabs1(dx) + Cabs1(e[0]) * Cabs1(x[j * ldx + 2]);
		for (i = 1; i < n - 1; i++) {
		    bi = B[i + j * ldb];
		    cx = conj(e[i - 1]) * x[i - 1 + j * ldx];
		    dx = d[i] * x[i + j * ldx];
		    ex = e[i] * x[i + 1 + j * ldx];
		    work[i] = bi - cx - dx - ex;
		    rwork[i] = Cabs1(bi) + Cabs1(e[i - 1]) * Cabs1(x[j * ldx + i - 1]) + Cabs1(dx) + Cabs1(e[i]) * Cabs1(x[j * ldx + i + 1]);
		}
		bi = B[n + j * ldb];
		cx = conj(e[n - 1]) * x[n - 1 + j * ldx];
		dx = d[n] * x[n + j * ldx];
		work[n] = bi - cx - dx;
		rwork[n] = Cabs1(bi) + Cabs1(e[n - 1]) * Cabs1(x[j * ldx + n - 1]) + Cabs1(dx);
	    }
	} else {
	    if (n == 1) {
		bi = B[j * ldb + 1];
		dx = d[1] * x[j * ldx + 1];
		work[1] = bi - dx;
		rwork[1] = Cabs1(bi) + Cabs1(dx);
	    } else {
		bi = B[j * ldb + 1];
		dx = d[1] * x[j * ldx + 1];
		ex = conj(e[0]) * x[j * ldx + 2];
		work[1] = bi - dx - ex;
		rwork[1] = Cabs1(bi) + Cabs1(dx) + Cabs1(e[0]) * Cabs1(x[j * ldx + 2]);
		for (i = 1; i < n - 1; i++) {
		    bi = B[i + j * ldb];
		    cx = e[i - 1] * x[i - 1 + j * ldx];
		    dx = d[i] * x[i + j * ldx];
		    ex = conj(e[i]) * x[i + 1 + j * ldx];
		    work[i] = bi - cx - dx - ex;
		    rwork[i] = Cabs1(bi) + Cabs1(e[i - 1]) * Cabs1(x[j * ldx + i - 1]) + Cabs1(dx) + Cabs1(e[i]) + Cabs1(x[j * ldx + i + 1]);
		}
		bi = B[n + j * ldb];
		cx = e[n - 1] * x[n - 1 + j * ldx];
		dx = d[n] * x[n + j * ldx];
		work[n] = bi - cx - dx;
		rwork[n] = Cabs1(bi) + Cabs1(e[n - 1]) * Cabs1(x[j * ldx + n - 1]) + Cabs1(dx);
	    }
	}
//Compute componentwise relative backward error from formula
//max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
//where abs(Z) is the componentwise absolute value of the matrix
//or vector Z.  If the i-th component of the denominator is less
//than SAFE2, then SAFE1 is added to the i-th components of the
//numerator and denominator before dividing.
	s = Zero;
	for (i = 0; i < n; i++) {
	    if (rwork[i] > safe2) {
		mtemp1 = s, mtemp2 = Cabs1(work[i]) / rwork[i];
		s = max(mtemp1, mtemp2);
	    } else {
		mtemp1 = s, mtemp2 = (Cabs1(work[i]) + safe1) / (rwork[i] + safe1);
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
	    Cpttrs(uplo, n, 1, &df[1], &ef[1], &work[0], n, info);
	    Caxpy(n, One, &work[0], 1, &x[j * ldx + 1], 1);
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
//Bound error from formula
//norm(X - XTRUE) / norm(X) .le. FERR =
//norm( abs(inv(A))*
//   ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)
// where
//  norm(Z) is the magnitude of the largest component of Z
//  inv(A) is the inverse of A
//  abs(Z) is the componentwise absolute value of the matrix or
//    vector Z
//  NZ is the maximum number of nonzeros in any row of A, plus 1
//  EPS is machine epsilon
//The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
//is incremented by SAFE1 if the i-th component of
//abs(A)*abs(X) + abs(B) is less than SAFE2
	for (i = 0; i < n; i++) {
	    if (rwork[i] > safe2) {
		rwork[i] = One + rwork[i - 1] * abs(ef[i - 1]);
	    } else {
		rwork[i] = rwork[i] / df[i] + rwork[i + 1] * abs(ef[i]);
	    }
	}
	ix = iRamax(n, &rwork[1], 1);
	ferr[j] = rwork[ix];
//Estimate the norm of inv(A).
//Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
//   m(i,j) =  abs(A(i,j)), i = j,
//   m(i,j) = -abs(A(i,j)), i .ne. j,
//and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.
//Solve M(L) * x = e.
	rwork[1] = One;
	for (i = 1; i < n; i++) {
	    rwork[i] = rwork[i - 1] * abs(ef[i - 1]) + One;
	}
//Solve D * M(L)' * x = b.
	rwork[n] = rwork[n] / df[n];
	for (i = n - 1; i >= 1; i--) {
	    rwork[i] = rwork[i] / df[i] + rwork[i + 1] * abs(ef[i]);
	}
//Compute norm(inv(A)) = max(x(i)), 1<=i<=n.
	ix = iRamax(n, &rwork[1], 1);
	ferr[j] = ferr[j] * abs(rwork[ix]);
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
