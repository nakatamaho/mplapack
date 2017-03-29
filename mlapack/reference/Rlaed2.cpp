/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed2.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaed2(INTEGER * k, INTEGER n, INTEGER n1, REAL * d, REAL * q, INTEGER ldq,
       INTEGER * indxq, REAL * rho, REAL * z, REAL * dlamda, REAL * w, REAL * q2, INTEGER * indx, INTEGER * indxc, INTEGER * indxp, INTEGER * coltyp, INTEGER * info)
{
    REAL c;
    INTEGER i, j;
    REAL s, t;
    INTEGER k2, n2, ct, nj, pj = 0, js, iq1, iq2, n1p1;
    REAL eps, tau, tol;
    INTEGER psm[4], imax, jmax;
    INTEGER ctot[4];
    REAL One = 1.0, Two = 2.0, Eight = 8.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    if (n < 0) {
	*info = -2;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -6;
    } else {
	if (max((INTEGER) 1, n / 2) > n1 || n / 2 < n1) {
	    *info = -3;
	}
    }
    if (*info != 0) {
	Mxerbla("Rlaed2", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

    n2 = n - n1;
    n1p1 = n1 + 1;

    if (*rho < Zero)
	Rscal(n2, -One, &z[n1p1], 1);

//Normalize z so that norm(z) = One  Since z is the concatenation of
//two normalized vectors, norm2(z) = sqrt(2).

    t = One / sqrt(Two);
    Rscal(n, t, &z[0], 1);

//RHO = ABS( norm(z)**2 * RHO )
    *rho = abs(*rho * Two);

//Sort the eigenvalues into increasing order
    for (i = n1p1; i < n; i++) {
	indxq[i] += n1;
    }
//re-integrate the deflated parts from the last pass
    for (i = 0; i < n; i++) {
	dlamda[i] = d[indxq[i]];
    }
    Rlamrg(n1, n2, &dlamda[0], 1, 1, &indxc[1]);
    for (i = 0; i < n; i++) {
	indx[i] = indxq[indxc[i]];

    }
//Calculate the allowable deflation tolerance
    imax = iRamax(n, &z[1], 1);
    jmax = iRamax(n, &d[0], 1);
    eps = Rlamch("Epsilon");
    mtemp1 = abs(d[jmax]), mtemp2 = abs(z[imax]);
    tol = eps * Eight * max(mtemp1, mtemp2);

//If the rank-1 modifier is small enough, no more needs to be done
//except to reorganize Q so that its columns correspond with the
//elements in D.

    if (*rho * abs(z[imax]) <= tol) {
	k = 0;
	iq2 = 1;
	for (j = 0; j < n; j++) {
	    i = indx[j];
	    Rcopy(n, &q[i * ldq + 1], 1, &q2[iq2], 1);
	    dlamda[j] = d[i];
	    iq2 += n;
	}
	Rlacpy("A", n, n, &q2[1], n, &q[0], ldq);
	Rcopy(n, &dlamda[1], 1, &d[0], 1);
	return;
    }
//If there are multiple eigenvalues then the problem deflates.  Here
//the number of equal eigenvalues are found.  As each equal
//eigenvalue is found, an elementary reflector is computed to rotate
//the corresponding eigensubspace so that the corresponding
//components of Z are zero in this new basis.

    for (i = 0; i < n; i++) {
	coltyp[i] = 1;

    }
    for (i = n1p1; i < n; i++) {
	coltyp[i] = 3;
    }

    k = 0;
    k2 = n + 1;
    for (j = 0; j < n; j++) {
	nj = indx[j];
	if (*rho * abs(z[nj]) <= tol) {
//Deflate due to small z component.
	    k2--;
	    coltyp[nj] = 4;
	    indxp[k2] = nj;
	    if (j == n) {
		goto L100;
	    }
	} else {
	    pj = nj;
	    break;
	}

    }
    j++;
    nj = indx[j];
    if (j > n) {
	goto L100;
    }
    if (*rho * abs(z[nj]) <= tol) {
// Deflate due to small z component.

	k2--;
	coltyp[nj] = 4;
	indxp[k2] = nj;
    } else {
//Check if eigenvalues are close enough to allow deflation.
	s = z[pj];
	c = z[nj];
//Find sqrt(a**2+b**2) without overflow or
//destructive underflow.
	tau = Rlapy2(c, s);
	t = d[nj] - d[pj];
	c /= tau;
	s = -s / tau;
	if (abs(t * c * s) <= tol) {
//Deflation is possible.
	    z[nj] = tau;
	    z[pj] = Zero;
	    if (coltyp[nj] != coltyp[pj]) {
		coltyp[nj] = 2;
	    }
	    coltyp[pj] = 4;
	    Rrot(n, &q[pj * ldq + 1], 1, &q[nj * ldq + 1], 1, c, s);
	    t = d[pj] * (c * c) + d[nj] * (s * s);
	    d[nj] = d[pj] * (s * s) + d[nj] * (c * c);
	    d[pj] = t;
	    k2--;
	    i = 1;
	    while (1) {
		if (k2 + i <= n) {
		    if (d[pj] < d[indxp[k2 + i]]) {
			indxp[k2 + i - 1] = indxp[k2 + i];
			indxp[k2 + i] = pj;
			i++;
			continue;
		    } else {
			indxp[k2 + i - 1] = pj;
		    }
		} else {
		    indxp[k2 + i - 1] = pj;
		}
	    }
	    pj = nj;
	} else {
	    ++(*k);
	    dlamda[*k] = d[pj];
	    w[*k] = z[pj];
	    indxp[*k] = pj;
	    pj = nj;
	}
    }
  L100:

//Record the last eigenvalue.
    ++(*k);
    dlamda[*k] = d[pj];
    w[*k] = z[pj];
    indxp[*k] = pj;
//Count up the total number of the various types of columns, then
//form a permutation which positions the four column types into
//four uniform groups (although one or more of these groups may be
//empty).
    for (j = 1; j < 4; j++) {
	ctot[j - 1] = 0;
    }
    for (j = 0; j < n; j++) {
	ct = coltyp[j];
	++ctot[ct - 1];
    }

//PSM(*) = Position in SubMatrix (of types 1 through 4)
    psm[0] = 1;
    psm[1] = ctot[0] + 1;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];
    *k = n - ctot[3];

//Fill out the INDXC array so that the permutation which it induces
//will place all type-1 columns first, all type-2 columns next,
//then all type-3's, and finally all type-4's.
    for (j = 0; j < n; j++) {
	js = indxp[j];
	ct = coltyp[js];
	indx[psm[ct - 1]] = js;
	indxc[psm[ct - 1]] = j;
	++psm[ct - 1];

    }

//Sort the eigenvalues and corresponding eigenvectors into DLAMDA
//and Q2 respectively.  The eigenvalues/vectors which were not
//deflated go into the first K slots of DLAMDA and Q2 respectively,
//while those which were deflated go into the last N - K slots.
    i = 1;
    iq1 = 1;
    iq2 = (ctot[0] + ctot[1]) * n1 + 1;
    for (j = 0; j < ctot[0]; j++) {
	js = indx[i];
	Rcopy(n1, &q[js * ldq + 1], 1, &q2[iq1], 1);
	z[i] = d[js];
	i++;
	iq1 += n1;

    }

    for (j = 0; j < ctot[1]; j++) {
	js = indx[i];
	Rcopy(n1, &q[js * ldq + 1], 1, &q2[iq1], 1);
	Rcopy(n2, &q[n1 + 1 + js * ldq], 1, &q2[iq2], 1);
	z[i] = d[js];
	i++;
	iq1 += n1;
	iq2 += n2;
    }

    for (j = 0; j < ctot[2]; j++) {
	js = indx[i];
	Rcopy(n2, &q[n1 + 1 + js * ldq], 1, &q2[iq2], 1);
	z[i] = d[js];
	i++;
	iq2 += n2;
    }
    iq1 = iq2;
    for (j = 0; j < ctot[3]; j++) {
	js = indx[i];
	Rcopy(n, &q[js * ldq + 1], 1, &q2[iq2], 1);
	iq2 += n;
	z[i] = d[js];
	i++;
    }

//The deflated eigenvalues and their corresponding vectors go back
//into the last N - K slots of D and Q respectively.
    Rlacpy("A", n, ctot[3], &q2[iq1], n, &q[(*k + 1) * ldq + 1], ldq);
    Rcopy(n - (*k), &z[*k + 1], 1, &d[*k + 1], 1);

//Copy CTOT into COLTYP for referencing in DLAED3.
    for (j = 0; j < 4; j++) {
	coltyp[j] = ctot[j - 1];
    }
    return;
}
