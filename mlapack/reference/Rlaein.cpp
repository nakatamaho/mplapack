/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaein.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaein(INTEGER rightv, INTEGER noinit, INTEGER n,
       REAL * h, INTEGER ldh, REAL wr, REAL wi, REAL * vr, REAL * vi, REAL * B, INTEGER ldb, REAL * work, REAL eps3, REAL smlnum, REAL bignum, INTEGER * info)
{
    INTEGER i, j;
    REAL w, x, y;
    INTEGER i1, i2, i3;
    REAL w1, ei, ej, xi, xr, rec;
    INTEGER its, ierr;
    REAL temp, norm, vmax;
    REAL scale;
    char trans;
    REAL vcrit, rootn, vnorm;
    REAL absbii, absbjj;
    char normin;
    REAL nrmsml, growto;
    REAL mtemp1, mtemp2;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
//GROWTO is the threshold used in the acceptance test for an
//eigenvector.
    rootn = (REAL) sqrt(n);
    growto = .1 / rootn;

    mtemp1 = One, mtemp2 = eps3 * rootn;
    nrmsml = max(mtemp1, mtemp2) * smlnum;

//Form B = H - (WR,WI)*I (except that the subdiagonal elements and
//the imaginary parts of the diagonal elements are not stored).
    for (j = 0; j < n; j++) {
	for (i = 0; i < j - 1; i++) {
	    B[i + j * ldb] = h[i + j * ldh];
	}
	B[j + j * ldb] = h[j + j * ldh] - wr;
    }
    if (wi == Zero) {
//Real eigenvalue.
	if (noinit) {
//Set initial vector.
	    for (i = 0; i < n; i++) {
		vr[i] = eps3;
	    }
	} else {
//Scale supplied initial vector.
	    vnorm = Rnrm2(n, &vr[1], 1);
	    mtemp1 = eps3 * rootn / max(vnorm, nrmsml);
	    Rscal(n, mtemp1, &vr[1], 1);
	}
	if (rightv) {
//LU decomposition with partial pivoting of B, replacing zero
//pivots by EPS3.
	    for (i = 0; i < n - 1; i++) {
		ei = h[i + 1 + i * ldh];
		if (abs(B[i + i * ldb]) < abs(ei)) {
//Interchange rows and eliminate.
		    x = B[i + i * ldb] / ei;
		    B[i + i * ldb] = ei;
		    for (j = i + 1; j <= n; j++) {
			temp = B[i + 1 + j * ldb];
			B[i + 1 + j * ldb] = B[i + j * ldb] - x * temp;
			B[i + j * ldb] = temp;
		    }
		} else {
//Eliminate without interchange.
		    if (B[i + i * ldb] == Zero) {
			B[i + i * ldb] = eps3;
		    }
		    x = ei / B[i + i * ldb];
		    if (x != Zero) {
			for (j = i + 1; j <= n; j++) {
			    B[i + 1 + j * ldb] = B[i + 1 + j * ldb] - x * B[i + j * ldb];
			}
		    }
		}
	    }
	    if (B[n + n * ldb] == Zero) {
		B[n + n * ldb] = eps3;
	    }
	    trans = 'N';
	} else {
//UL decomposition with partial pivoting of B, replacing zero
//pivots by EPS3.
	    for (j = n; j >= 2; j--) {
		ej = h[j + (j - 1) * ldh];
		if (abs(B[j + j * ldb]) < abs(ej)) {
//Interchange columns and eliminate.
		    x = B[j + j * ldb] / ej;
		    B[j + j * ldb] = ej;
		    for (i = 0; i < j - 1; i++) {
			temp = B[i + (j - 1) * ldb];
			B[i + (j - 1) * ldb] = B[i + j * ldb] - x * temp;
			B[i + j * ldb] = temp;
		    }
		} else {
//Eliminate without interchange.
		    if (B[j + j * ldb] == Zero) {
			B[j + j * ldb] = eps3;
		    }
		    x = ej / B[j + j * ldb];
		    if (x != Zero) {
			for (i = 0; i < j - 1; i++) {
			    B[i + (j - 1) * ldb] = B[i + (j - 1) * ldb] - x * B[i + j * ldb];
			}
		    }
		}
	    }
	    if (B[ldb + 1] == Zero) {
		B[ldb + 1] = eps3;
	    }
	    trans = 'T';
	}
	normin = 'N';
	for (its = 1; its <= n; its++) {
//Solve U*x = scale*v for a right eigenvector
//  or U'*x = scale*v for a left eigenvector,
//overwriting x on v.
	    Rlatrs("Upper", &trans, "Nonunit", &normin, n, &B[0], ldb, &vr[1], &scale, &work[0], &ierr);
	    normin = 'Y';
//Test for sufficient growth in the norm of v.
	    vnorm = Rasum(n, &vr[1], 1);
	    if (vnorm >= growto * scale) {
		goto L120;
	    }
//Choose new orthogonal starting vector and try again.
	    temp = eps3 / (rootn + One);
	    vr[1] = eps3;
	    for (i = 1; i < n; i++) {
		vr[i] = temp;
	    }
	    vr[n - its + 1] -= eps3 * rootn;
	}
//Failure to find eigenvector in N iterations.
	*info = 1;
      L120:
//Normalize eigenvector.
	i = iRamax(n, &vr[1], 1);
	Rscal(n, One / abs(vr[i]), &vr[1], 1);
    } else {
//Complex eigenvalue.
	if (noinit) {
//Set initial vector.
	    for (i = 0; i < n; i++) {
		vr[i] = eps3;
		vi[i] = Zero;
	    }
	} else {
//Scale supplied initial vector.
	    mtemp1 = Rnrm2(n, &vr[1], 1);
	    mtemp2 = Rnrm2(n, &vi[1], 1);
	    norm = Rlapy2(mtemp1, mtemp2);
	    rec = eps3 * rootn / max(norm, nrmsml);
	    Rscal(n, rec, &vr[1], 1);
	    Rscal(n, rec, &vi[1], 1);
	}

	if (rightv) {
//LU decomposition with partial pivoting of B, replacing zero
//pivots by EPS3.
//The imaginary part of the (i,j)-th element of U is stored in
//B(j+1,i).
	    B[ldb + 2] = -(wi);
	    for (i = 1; i < n; i++) {
		B[i + 1 + ldb] = Zero;
	    }
	    for (i = 0; i < n - 1; i++) {
		absbii = Rlapy2(B[i + i * ldb], B[i + 1 + i * ldb]);
		ei = h[i + 1 + i * ldh];
		if (absbii < abs(ei)) {
//Interchange rows and eliminate.
		    xr = B[i + i * ldb] / ei;
		    xi = B[i + 1 + i * ldb] / ei;
		    B[i + i * ldb] = ei;
		    B[i + 1 + i * ldb] = Zero;
		    for (j = i + 1; j <= n; j++) {
			temp = B[i + 1 + j * ldb];
			B[i + 1 + j * ldb] = B[i + j * ldb] - xr * temp;
			B[j + 1 + (i + 1) * ldb] = B[j + 1 + i * ldb] - xi * temp;
			B[i + j * ldb] = temp;
			B[j + 1 + i * ldb] = Zero;
		    }
		    B[i + 2 + i * ldb] = -(wi);
		    B[i + 1 + (i + 1) * ldb] -= xi * wi;
		    B[i + 2 + (i + 1) * ldb] += xr * wi;
		} else {
//Eliminate without interchanging rows.
		    if (absbii == Zero) {
			B[i + i * ldb] = eps3;
			B[i + 1 + i * ldb] = Zero;
			absbii = eps3;
		    }
		    ei = ei / absbii / absbii;
		    xr = B[i + i * ldb] * ei;
		    xi = -B[i + 1 + i * ldb] * ei;
		    for (j = i + 1; j <= n; j++) {
			B[i + 1 + j * ldb] = B[i + 1 + j * ldb] - xr * B[i + j * ldb] + xi * B[j + 1 + i * ldb];
			B[j + 1 + (i + 1) * ldb] = -xr * B[j + 1 + i * ldb] - xi * B[i + j * ldb];
		    }
		    B[i + 2 + (i + 1) * ldb] -= wi;
		}
//Compute 1-norm of offdiagonal elements of i-th row.
		work[i] = Rasum(n - i, &B[i + (i + 1) * ldb], ldb)
		    + Rasum(n - i, &B[i + 2 + i * ldb], 1);
	    }
	    if (B[n + n * ldb] == Zero && B[n + 1 + n * ldb] == Zero) {
		B[n + n * ldb] = eps3;
	    }
	    work[n] = Zero;
	    i1 = n;
	    i2 = 1;
	    i3 = -1;
	} else {
//UL decomposition with partial pivoting of conjg(B),
//replacing zero pivots by EPS3.
//The imaginary part of the (i,j)-th element of U is stored in
//B(j+1,i).
	    B[n + 1 + n * ldb] = wi;
	    for (j = 0; j < n - 1; j++) {
		B[n + 1 + j * ldb] = Zero;
	    }
	    for (j = n; j >= 2; j--) {
		ej = h[j + (j - 1) * ldh];
		absbjj = Rlapy2(B[j + j * ldb], B[j + 1 + j * ldb]);
		if (absbjj < abs(ej)) {
//Interchange columns and eliminate
		    xr = B[j + j * ldb] / ej;
		    xi = B[j + 1 + j * ldb] / ej;
		    B[j + j * ldb] = ej;
		    B[j + 1 + j * ldb] = Zero;
		    for (i = 0; i < j - 1; i++) {
			temp = B[i + (j - 1) * ldb];
			B[i + (j - 1) * ldb] = B[i + j * ldb] - xr * temp;
			B[j + i * ldb] = B[j + 1 + i * ldb] - xi * temp;
			B[i + j * ldb] = temp;
			B[j + 1 + i * ldb] = Zero;
		    }
		    B[j + 1 + (j - 1) * ldb] = wi;
		    B[j - 1 + (j - 1) * ldb] = B[j - 1 + (j - 1) * ldb] + xi * wi;
		    B[j + (j - 1) * ldb] = B[j + (j - 1) * ldb] - xr * wi;
		} else {
//Eliminate without interchange.
		    if (absbjj == Zero) {
			B[j + j * ldb] = eps3;
			B[j + 1 + j * ldb] = Zero;
			absbjj = eps3;
		    }
		    ej = ej / absbjj / absbjj;
		    xr = B[j + j * ldb] * ej;
		    xi = -B[j + 1 + j * ldb] * ej;
		    for (i = 0; i < j - 1; i++) {
			B[i + (j - 1) * ldb] = B[i + (j - 1) * ldb]
			    - xr * B[i + j * ldb] + xi * B[j + 1 + i * ldb];
			B[j + i * ldb] = -xr * B[j + 1 + i * ldb] - xi * B[i + j * ldb];

		    }
		    B[j + (j - 1) * ldb] += wi;
		}
//Compute 1-norm of offdiagonal elements of j-th column.
		work[j] = Rasum(j - 1, &B[j * ldb + 1], 1) + Rasum(j - 1, &B[j + 1 + ldb], ldb);
	    }
	    if (B[ldb + 1] == Zero && B[ldb + 2] == Zero) {
		B[ldb + 1] = eps3;
	    }
	    work[1] = Zero;
	    i1 = 1;
	    i2 = n;
	    i3 = 1;
	}
	for (its = 1; its <= n; its++) {
	    scale = One;
	    vmax = One;
	    vcrit = bignum;
//Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
//  or U'*(xr,xi) = scale*(vr,vi) for a left eigenvector,
//overwriting (xr,xi) on (vr,vi).
	    for (i = i1; i <= i2; i += i3) {
		if (work[i] > vcrit) {
		    rec = One / vmax;
		    Rscal(n, rec, &vr[1], 1);
		    Rscal(n, rec, &vi[1], 1);
		    scale *= rec;
		    vmax = One;
		    vcrit = bignum;
		}
		xr = vr[i];
		xi = vi[i];
		if (rightv) {
		    for (j = i + 1; j <= n; j++) {
			xr = xr - B[i + j * ldb] * vr[j] + B[j + 1 + i * ldb] * vi[j];
			xi = xi - B[i + j * ldb] * vi[j] - B[j + 1 + i * ldb] * vr[j];

		    }
		} else {
		    for (j = 0; j < i - 1; j++) {
			xr = xr - B[j + i * ldb] * vr[j] + B[i + 1 + j * ldb] * vi[j];
			xi = xi - B[j + i * ldb] * vi[j] - B[i + 1 + j * ldb] * vr[j];

		    }
		}
		w = abs(B[i + i * ldb]) + abs(B[i + 1 + i * ldb]);
		if (w > smlnum) {
		    if (w < One) {
			w1 = abs(xr) + abs(xi);
			if (w1 > w * bignum) {
			    rec = One / w1;
			    Rscal(n, rec, &vr[1], 1);
			    Rscal(n, rec, &vi[1], 1);
			    xr = vr[i];
			    xi = vi[i];
			    scale *= rec;
			    vmax *= rec;
			}
		    }
//Divide by diagonal element of B.
		    Rladiv(xr, xi, B[i + i * ldb], B[i + 1 + i * ldb], &vr[i], &vi[i]);
		    mtemp1 = abs(vr[i]) + abs(vi[i]);
		    mtemp2 = vmax;
		    vmax = max(mtemp1, mtemp2);
		    vcrit = bignum / vmax;
		} else {
		    for (j = 0; j < n; j++) {
			vr[j] = Zero;
			vi[j] = Zero;

		    }
		    vr[i] = One;
		    vi[i] = One;
		    scale = Zero;
		    vmax = One;
		    vcrit = bignum;
		}
	    }
//Test for sufficient growth in the norm of (VR,VI).
	    vnorm = Rasum(n, &vr[1], 1) + Rasum(n, &vi[1], 1);
	    if (vnorm >= growto * scale) {
		goto L280;
	    }
//Choose a new orthogonal starting vector and try again.
	    y = eps3 / (rootn + One);
	    vr[1] = eps3;
	    vi[1] = Zero;
	    for (i = 1; i < n; i++) {
		vr[i] = y;
		vi[i] = Zero;
	    }
	    vr[n - its + 1] = vr[n - its + 1] - eps3 * rootn;
	}
//Failure to find eigenvector in N iterations
	*info = 1;
      L280:
//Normalize eigenvector.
	vnorm = Zero;
	for (i = 0; i < n; i++) {
	    mtemp1 = vnorm, mtemp2 = abs(vr[i]) + abs(vi[i]);
	    vnorm = max(mtemp1, mtemp2);
	}
	Rscal(n, One / vnorm, &vr[1], 1);
	Rscal(n, One / vnorm, &vi[1], 1);

    }
    return;
}
