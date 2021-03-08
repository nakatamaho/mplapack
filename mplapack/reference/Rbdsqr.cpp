/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rbdsqr.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rbdsqr(const char *uplo, INTEGER n, INTEGER ncvt, INTEGER nru, INTEGER ncc, REAL * d,
       REAL * e, REAL * vt, INTEGER ldvt, REAL * u, INTEGER ldu, REAL * c, INTEGER ldc, REAL * work, INTEGER * info)
{

    REAL f, g, h;
    INTEGER i, j, m;
    REAL r, cs;
    INTEGER ll;
    REAL sn, mu;
    INTEGER nm1, nm12, nm13, lll;
    REAL eps, sll, tol, abse;
    INTEGER idir;
    REAL abss;
    INTEGER oldm;
    REAL cosl;
    INTEGER isub, iter;
    REAL unfl, sinl, cosr, smin, smax, sinr;
    REAL oldcs;
    INTEGER oldll;
    REAL shift, sigmn, oldsn;
    INTEGER maxit;
    REAL sminl, sigmx;
    INTEGER lower;
    REAL sminoa, thresh;
    INTEGER rotate;
    REAL tolmul;
    REAL mtemp1, mtemp2;
    REAL OneEighth = 0.125;
    REAL One = 1.0, Zero = 0.0, Ten = 10.0, Hndrd = 100.0;


//Test the input parameters.
    *info = 0;
    lower = Mlsame(uplo, "L");
    if (!Mlsame(uplo, "U") && !lower) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (ncvt < 0) {
	*info = -3;
    } else if (nru < 0) {
	*info = -4;
    } else if (ncc < 0) {
	*info = -5;
    } else if ((ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < max((INTEGER) 1, n))) {
	*info = -9;
    } else if (ldu < max((INTEGER) 1, nru)) {
	*info = -11;
    } else if ((ncc == 0 && ldc < 1) || (ncc > 0 && ldc < max((INTEGER) 1, n))) {
	*info = -13;
    }
    if (*info != 0) {
	Mxerbla("Rdbsqr", -(*info));
	return;
    }
    if (n == 0) {
	return;
    }
    if (n == 1) {
	goto L160;
    }
//ROTATE is true if any singular vectors desired, false otherwise
    rotate = ncvt > 0 || nru > 0 || ncc > 0;

//If no singular vectors desired, use qd algorithm
    if (!rotate) {
	Rlasq1(n, &d[0], &e[0], &work[0], info);
	return;
    }
    nm1 = n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;

//Get machine constants
    eps = Rlamch("Epsilon");
    unfl = Rlamch("Safe minimum");

//If matrix lower bidiagonal, rotate to be upper bidiagonal
//by applying Givens rotations on the left
    if (lower) {
	for (i = 0; i < n - 1; i++) {
	    Rlartg(d[i], e[i], &cs, &sn, &r);
	    d[i] = r;
	    e[i] = sn * d[i + 1];
	    d[i + 1] = cs * d[i + 1];
	    work[i] = cs;
	    work[nm1 + i] = sn;
	}
//Update singular vectors if desired
	if (nru > 0) {
	    Rlasr("R", "V", "F", nru, n, &work[0], &work[n], &u[0], ldu);
	}
	if (ncc > 0) {
	    Rlasr("L", "V", "F", n, ncc, &work[0], &work[n], &c[0], ldc);
	}
    }
//Compute singular values to relative accuracy TOL
//(By setting TOL to be negative, algorithm will compute
//singular values to absolute accuracy ABS(TOL)norm(input matrix))
    mtemp1 = pow(eps, -OneEighth);
    mtemp2 = min(Hndrd, mtemp1);
    tolmul = max(Ten, mtemp2);
    tol = tolmul * eps;

//Compute approximate maximum, minimum singular values

    smax = Zero;
    for (i = 0; i < n; i++) {
	mtemp1 = smax, mtemp2 = abs(d[i]);
	smax = max(mtemp1, mtemp2);
    }
    for (i = 0; i < n - 1; i++) {
	mtemp1 = smax, mtemp2 = abs(e[i]);
	smax = max(mtemp1, mtemp2);
    }
    sminl = Zero;
    if (tol >= Zero) {
//Relative accuracy desired
	sminoa = abs(d[1]);
	if (sminoa != Zero) {
	    mu = sminoa;
	    for (i = 1; i < n; i++) {
		mu = abs(d[i]) * (mu / (mu + abs(e[i - 1])));
		sminoa = min(sminoa, mu);
		if (sminoa == Zero)
		    break;
	    }
	}
	sminoa /= sqrt(n);
	mtemp1 = tol * sminoa, mtemp2 = n * 6 * n * unfl;
	thresh = max(mtemp1, mtemp2);
    } else {
//Absolute accuracy desired
	mtemp1 = abs(tol) * smax;
	mtemp2 = n * 6 * n * unfl;
	thresh = max(mtemp1, mtemp2);
    }

//Prepare for main iteration loop for the singular values
//(MAXIT is the maximum number of passes through the inner
//loop permitted before nonconvergence signalled.)
    maxit = n * 6 * n;
    iter = 0;
    oldll = -1;
    oldm = -1;
//M points to last element of unconverged part of matrix
    m = n;
//Begin main iteration loop
  L60:
//Check for convergence or exceeding iteration count
    if (m <= 1) {
	goto L160;
    }
    if (iter > maxit) {
	goto L200;
    }
//Find diagonal block of matrix to work on
    if (tol < Zero && abs(d[m]) <= thresh) {
	d[m] = Zero;
    }
    smax = abs(d[m]);
    smin = smax;
    for (lll = 0; lll <= m - 1; lll++) {
	ll = m - lll;
	abss = abs(d[ll]);
	abse = abs(e[ll]);
	if (tol < Zero && abss <= thresh) {
	    d[ll] = Zero;
	}
	if (abse <= thresh) {
	    goto L80;
	}
	smin = min(smin, abss);
	mtemp1 = max(smax, abss);
	smax = max(mtemp1, abse);
    }
    ll = 0;
    goto L90;
  L80:
    e[ll] = Zero;
//Matrix splits since E(LL) = 0
    if (ll == m - 1) {
//Convergence of bottom singular value, return to top of loop
	m--;
	goto L60;
    }
  L90:
    ll++;
//E(LL) through E(M-1) are nonzero, E(LL-1) is zero
    if (ll == m - 1) {
//2 by 2 block, handle separately
	Rlasv2(d[m - 1], e[m - 1], d[m], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
	d[m - 1] = sigmx;
	e[m - 1] = Zero;
	d[m] = sigmn;
//Compute singular vectors, if desired
	if (ncvt > 0) {
	    Rrot(ncvt, &vt[m - 1 + ldvt], ldvt, &vt[m + ldvt], ldvt, cosr, sinr);
	}
	if (nru > 0) {
	    Rrot(nru, &u[(m - 1) * ldu + 1], 1, &u[m * ldu + 1], 1, cosl, sinl);
	}
	if (ncc > 0) {
	    Rrot(ncc, &c[m - 1 + ldc], ldc, &c[m + ldc], ldc, cosl, sinl);
	}
	m += -2;
	goto L60;
    }
//If working on new submatrix, choose shift direction
//(from larger end diagonal element towards smaller)
    if (ll > oldm || m < oldll) {
	if (abs(d[ll]) >= abs(d[m])) {
//Chase bulge from top (big end) to bottom (small end)
	    idir = 1;
	} else {
//Chase bulge from bottom (big end) to top (small end)
	    idir = 2;
	}
    }
// Apply convergence tests
    if (idir == 1) {
//Run convergence test in forward direction
//First apply standard test to bottom of matrix
	if ((abs(e[m - 1]) <= abs(tol) * abs(d[m])) || (tol < Zero && abs(e[m - 1]) <= thresh)) {
	    e[m - 1] = Zero;
	    goto L60;
	}
	if (tol >= Zero) {
//If relative accuracy desired,
//apply convergence criterion forward
	    mu = abs(d[ll]);
	    sminl = mu;
	    for (lll = ll; lll <= m - 1; lll++) {
		if (abs(e[lll]) <= tol * mu) {
		    e[lll] = Zero;
		    goto L60;
		}
		mu = abs(d[lll + 1]) * (mu / (mu + abs(e[lll])));
		sminl = min(sminl, mu);
	    }
	}
    } else {
//Run convergence test in backward direction
//First apply standard test to top of matrix

	if ((abs(e[ll]) <= abs(tol) * abs(d[ll])) || (tol < Zero && abs(e[ll]) <= thresh)) {
	    e[ll] = Zero;
	    goto L60;
	}
	if (tol >= Zero) {
//If relative accuracy desired,
//apply convergence criterion backward

	    mu = abs(d[m]);
	    sminl = mu;
	    for (lll = m - 1; lll >= ll; lll--) {
		if (abs(e[lll]) <= tol * mu) {
		    e[lll] = Zero;
		    goto L60;
		}
		mu = abs(d[lll]) * (mu / (mu + abs(e[lll])));
		sminl = min(sminl, mu);
	    }
	}
    }
    oldll = ll;
    oldm = m;
//Compute shift.  First, test if shifting would ruin relative
//accuracy, and if so set the shift to zero.
    mtemp1 = eps, mtemp2 = tol * .01;
    if (tol >= Zero && n * tol * (sminl / smax) <= max(mtemp1, mtemp2)) {
//Use a zero shift to avoid loss of relative accuracy
	shift = Zero;
    } else {
//Compute the shift from 2-by-2 block at end of matrix
	if (idir == 1) {
	    sll = abs(d[ll]);
	    Rlas2(d[m - 1], e[m - 1], d[m], &shift, &r);
	} else {
	    sll = abs(d[m]);
	    Rlas2(d[ll], e[ll], d[ll + 1], &shift, &r);
	}
//Test if shift negligible, and if so set to zero
	if (sll > Zero) {
	    if ((shift / sll) * (shift / sll) < eps) {
		shift = Zero;
	    }
	}
    }
//Increment iteration count
    iter = iter + m - ll;
//If SHIFT = 0, do simplified QR iteration
    if (shift == Zero) {
	if (idir == 1) {
//Chase bulge from top to bottom
//Save cosines and sines for later singular vector updates
	    cs = One;
	    oldcs = One;
	    for (i = ll; i <= m - 1; i++) {
		Rlartg(d[i] * cs, e[i], &cs, &sn, &r);
		if (i > ll) {
		    e[i - 1] = oldsn * r;
		}
		Rlartg(oldcs * r, d[i + 1] * sn, &oldcs, &oldsn, &d[i]);
		work[i - ll + 1] = cs;
		work[i - ll + 1 + nm1] = sn;
		work[i - ll + 1 + nm12] = oldcs;
		work[i - ll + 1 + nm13] = oldsn;

	    }
	    h = d[m] * cs;
	    d[m] = h * oldcs;
	    e[m - 1] = h * oldsn;
//Update singular vectors
	    if (ncvt > 0) {
		Rlasr("L", "V", "F", m - ll + 1, ncvt, &work[0], &work[n], &vt[ll + ldvt], ldvt);
	    }
	    if (nru > 0) {
		Rlasr("R", "V", "F", nru, m - ll + 1, &work[nm12 + 1], &work[nm13 + 1], &u[ll * ldu + 1], ldu);
	    }
	    if (ncc > 0) {
		Rlasr("L", "V", "F", m - ll + 1, ncc, &work[nm12 + 1], &work[nm13 + 1], &c[ll + ldc], ldc);
	    }
//Test convergence
	    if (abs(e[m - 1]) <= thresh) {
		e[m - 1] = Zero;
	    }
	} else {

//Chase bulge from bottom to top
//Save cosines and sines for later singular vector updates
	    cs = One;
	    oldcs = One;
	    for (i = m; i >= ll + 1; i--) {
		Rlartg(d[i] * cs, e[i - 1], &cs, &sn, &r);
		if (i < m) {
		    e[i] = oldsn * r;
		}
		Rlartg(oldcs * r, d[i - 1] * sn, &oldcs, &oldsn, &d[i]);
		work[i - ll] = cs;
		work[i - ll + nm1] = -sn;
		work[i - ll + nm12] = oldcs;
		work[i - ll + nm13] = -oldsn;
	    }
	    h = d[ll] * cs;
	    d[ll] = h * oldcs;
	    e[ll] = h * oldsn;
//Update singular vectors
	    if (ncvt > 0) {
		Rlasr("L", "V", "B", m - ll + 1, ncvt, &work[nm12 + 1], &work[nm13 + 1], &vt[ll + ldvt], ldvt);
	    }
	    if (nru > 0) {
		Rlasr("R", "V", "B", nru, m - ll + 1, &work[0], &work[n], &u[ll * ldu + 1], ldu);
	    }
	    if (ncc > 0) {
		Rlasr("L", "V", "B", m - ll + 1, ncc, &work[0], &work[n], &c[ll + ldc], ldc);
	    }
//Test convergence
	    if (abs(e[ll]) <= thresh) {
		e[ll] = Zero;
	    }
	}
    } else {
//Use nonzero shift
	if (idir == 1) {
//Chase bulge from top to bottom
//Save cosines and sines for later singular vector updates
	    f = (abs(d[ll]) - shift) * (sign(One, d[ll]) + shift / d[ll]);
	    g = e[ll];
	    for (i = ll; i <= m - 1; i++) {
		Rlartg(f, g, &cosr, &sinr, &r);
		if (i > ll) {
		    e[i - 1] = r;
		}
		f = cosr * d[i] + sinr * e[i];
		e[i] = cosr * e[i] - sinr * d[i];
		g = sinr * d[i + 1];
		d[i + 1] = cosr * d[i + 1];
		Rlartg(f, g, &cosl, &sinl, &r);
		d[i] = r;
		f = cosl * e[i] + sinl * d[i + 1];
		d[i + 1] = cosl * d[i + 1] - sinl * e[i];
		if (i < m - 1) {
		    g = sinl * e[i + 1];
		    e[i + 1] = cosl * e[i + 1];
		}
		work[i - ll + 1] = cosr;
		work[i - ll + 1 + nm1] = sinr;
		work[i - ll + 1 + nm12] = cosl;
		work[i - ll + 1 + nm13] = sinl;

	    }
	    e[m - 1] = f;
//Update singular vectors
	    if (ncvt > 0) {
		Rlasr("L", "V", "F", m - ll + 1, ncvt, &work[0], &work[n], &vt[ll + ldvt], ldvt);
	    }
	    if (nru > 0) {
		Rlasr("R", "V", "F", nru, m - ll + 1, &work[nm12 + 1], &work[nm13 + 1], &u[ll * ldu + 1], ldu);
	    }
	    if (ncc > 0) {
		Rlasr("L", "V", "F", m - ll + 1, ncc, &work[nm12 + 1], &work[nm13 + 1], &c[ll + ldc], ldc);
	    }
//Test convergence
	    if (abs(e[m - 1]) <= thresh) {
		e[m - 1] = Zero;
	    }
	} else {
//Chase bulge from bottom to top
//Save cosines and sines for later singular vector updates
	    f = (abs(d[m]) - shift) * (sign(One, d[m]) + shift / d[m]);
	    g = e[m - 1];
	    for (i = m; i >= ll + 1; i--) {
		Rlartg(f, g, &cosr, &sinr, &r);
		if (i < m) {
		    e[i] = r;
		}
		f = cosr * d[i] + sinr * e[i - 1];
		e[i - 1] = cosr * e[i - 1] - sinr * d[i];
		g = sinr * d[i - 1];
		d[i - 1] = cosr * d[i - 1];
		Rlartg(f, g, &cosl, &sinl, &r);
		d[i] = r;
		f = cosl * e[i - 1] + sinl * d[i - 1];
		d[i - 1] = cosl * d[i - 1] - sinl * e[i - 1];
		if (i > ll + 1) {
		    g = sinl * e[i - 2];
		    e[i - 2] = cosl * e[i - 2];
		}
		work[i - ll] = cosr;
		work[i - ll + nm1] = -sinr;
		work[i - ll + nm12] = cosl;
		work[i - ll + nm13] = -sinl;

	    }
	    e[ll] = f;
//Test convergence
	    if (abs(e[ll]) <= thresh) {
		e[ll] = Zero;
	    }
//Update singular vectors if desired
	    if (ncvt > 0) {
		Rlasr("L", "V", "B", m - ll + 1, ncvt, &work[nm12 + 1], &work[nm13 + 1], &vt[ll + ldvt], ldvt);
	    }
	    if (nru > 0) {
		Rlasr("R", "V", "B", nru, m - ll + 1, &work[0], &work[n], &u[ll * ldu + 1], ldu);
	    }
	    if (ncc > 0) {
		Rlasr("L", "V", "B", m - ll + 1, ncc, &work[0], &work[n], &c[ll + ldc], ldc);
	    }
	}
    }

//QR iteration finished, go back and check convergence
    goto L60;

//All singular values converged, so make them positive
  L160:
    for (i = 0; i < n; i++) {
	if (d[i] < Zero) {
	    d[i] = -d[i];
//Change sign of singular vectors, if desired
	    if (ncvt > 0) {
		Rscal(ncvt, -One, &vt[i + ldvt], ldvt);
	    }
	}
    }
//Sort the singular values into decreasing order (insertion sort on
//singular values, but only one transposition per singular vector)
    for (i = 0; i < n - 1; i++) {
//Scan for smallest D(I)
	isub = 1;
	smin = d[1];
	for (j = 2; j <= n + 1 - i; j++) {
	    if (d[j] <= smin) {
		isub = j;
		smin = d[j];
	    }
	}
	if (isub != n + 1 - i) {
//Swap singular values and vectors
	    d[isub] = d[n + 1 - i];
	    d[n + 1 - i] = smin;
	    if (ncvt > 0) {
		Rswap(ncvt, &vt[isub + ldvt], ldvt, &vt[n + 1 - i + ldvt], ldvt);
	    }
	    if (nru > 0) {
		Rswap(nru, &u[isub * ldu + 1], 1, &u[(n + 1 - i) * ldu + 1], 1);
	    }
	    if (ncc > 0) {
		Rswap(ncc, &c[isub + ldc], ldc, &c[n + 1 - i + ldc], ldc);
	    }
	}
    }
    return;

//Maximum number of iterations exceeded, failure to converge
  L200:
    *info = 0;
    for (i = 0; i < n - 1; i++) {
	if (e[i] != Zero) {
	    ++(*info);
	}
    }
    return;
}
