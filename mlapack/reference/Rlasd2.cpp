/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd2.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasd2(INTEGER nl, INTEGER nr, INTEGER sqre, INTEGER * k, REAL * d, REAL * z,
       REAL alpha, REAL beta, REAL * u, INTEGER ldu, REAL * vt,
       INTEGER ldvt, REAL * dsigma, REAL * u2, INTEGER ldu2, REAL * vt2,
       INTEGER ldvt2, INTEGER * idxp, INTEGER * idx, INTEGER * idxc, INTEGER * idxq, INTEGER * coltyp, INTEGER * info)
{
    REAL c = 0.0;
    INTEGER i, j, m, n;
    REAL s = 0.0;
    INTEGER k2;
    REAL z1;
    INTEGER ct, jp;
    REAL eps, tau, tol;
    INTEGER psm[4], nlp1, nlp2, idxi, idxj;
    INTEGER ctot[4], idxjp;
    INTEGER jprev = 0;
    REAL hlftol;
    REAL Eight = 8.0, Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    if (nl < 1) {
	*info = -1;
    } else if (nr < 1) {
	*info = -2;
    } else if (sqre != 1 && sqre != 0) {
	*info = -3;
    }
    n = nl + nr + 1;
    m = n + sqre;
    if (ldu < n) {
	*info = -10;
    } else if (ldvt < m) {
	*info = -12;
    } else if (ldu2 < n) {
	*info = -15;
    } else if (ldvt2 < m) {
	*info = -17;
    }
    if (*info != 0) {
	Mxerbla("Rlasd2", -(*info));
	return;
    }
    nlp1 = nl + 1;
    nlp2 = nl + 2;
//Generate the first part of the vector Z; and move the singular
//values in the first part of D one position backward.
    z1 = alpha * vt[nlp1 + nlp1 * ldvt];
    z[1] = z1;
    for (i = nl; i >= 1; i--) {
	z[i + 1] = alpha * vt[i + nlp1 * ldvt];
	d[i + 1] = d[i];
	idxq[i + 1] = idxq[i] + 1;
    }
//Generate the second part of the vector Z.
    for (i = nlp2; i <= m; i++) {
	z[i] = beta * vt[i + nlp2 * ldvt];
    }
//Initialize some reference arrays.
    for (i = 1; i < nlp1; i++) {
	coltyp[i] = 1;
    }
    for (i = nlp2; i <= n; i++) {
	coltyp[i] = 2;
    }
//Sort the singular values into increasing order
    for (i = nlp2; i <= n; i++) {
	idxq[i] += nlp1;
    }
//DSIGMA, IDXC, IDXC, and the first column of U2
//are used as storage space.
    for (i = 1; i < n; i++) {
	dsigma[i] = d[idxq[i]];
	u2[i + ldu2] = z[idxq[i]];
	idxc[i] = coltyp[idxq[i]];
    }
    Rlamrg(nl, nr, &dsigma[2], 1, 1, &idx[2]);
    for (i = 1; i < n; i++) {
	idxi = idx[i] + 1;
	d[i] = dsigma[idxi];
	z[i] = u2[idxi + ldu2];
	coltyp[i] = idxc[idxi];
    }
//Calculate the allowable deflation tolerance
    eps = Rlamch("Epsilon");
    mtemp1 = max(abs(alpha), abs(beta));
    mtemp2 = abs(d[n]);
    tol = eps * Eight * max(mtemp1, mtemp2);
//There are 2 kinds of deflation -- first a value in the z-vector 
//is small, second two (or more) singular values are very close 
//together (their difference is small). 
//If the value in the z-vector is small, we simply permute the 
//array so that the corresponding singular value is moved to the 
//end. 
//If two values in the D-vector are close, we perform a two-sided 
//rotation designed to make one of the corresponding z-vector 
//entries zero, and then permute the array so that the deflated 
//singular value is moved to the end. 
//If there are multiple singular values then the problem deflates. 
//Here the number of equal singular values are found.  As each equal 
//singular value is found, an elementary reflector is computed to 
//rotate the corresponding singular subspace so that the 
//corresponding components of Z are zero in this new basis. 
    k = 0;
    k2 = n + 1;
    for (j = 2; j <= n; j++) {
	if (abs(z[j]) <= tol) {
//Deflate due to small z component.
	    k2--;
	    idxp[k2] = j;
	    coltyp[j] = 4;
	    if (j == n) {
		goto L120;
	    }
	} else {
	    jprev = j;
	    goto L90;
	}

    }
  L90:
    j = jprev;
  L100:
    j++;
    if (j > n) {
	goto L110;
    }
    if (abs(z[j]) <= tol) {
//Deflate due to small z component.
	k2--;
	idxp[k2] = j;
	coltyp[j] = 4;
    } else {
//Check if singular values are close enough to allow deflation.
	if (abs(d[j] - d[jprev]) <= tol) {
//Deflation is possible.
	    s = z[jprev];
	    c = z[j];
//Find sqrt(a**2+b**2) without overflow or
//destructive underflow.
	    tau = Rlapy2(c, s);
	    c /= tau;
	    s = -s / tau;
	    z[j] = tau;
	    z[jprev] = Zero;
//Apply back the Givens rotation to the left and right
//singular vector matrices.
	    idxjp = idxq[idx[jprev] + 1];
	    idxj = idxq[idx[j] + 1];
	    if (idxjp <= nlp1) {
		idxjp--;
	    }
	    if (idxj <= nlp1) {
		idxj--;
	    }
	    Rrot(n, &u[idxjp * ldu + 1], 1, &u[idxj * ldu + 1], 1, c, s);
	    Rrot(m, &vt[idxjp + ldvt], ldvt, &vt[idxj + ldvt], ldvt, c, s);
	    if (coltyp[j] != coltyp[jprev]) {
		coltyp[j] = 3;
	    }
	    coltyp[jprev] = 4;
	    k2--;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    u2[*k + ldu2] = z[jprev];
	    dsigma[*k] = d[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L100;
  L110:

//Record the last singular value.

    ++(*k);
    u2[*k + ldu2] = z[jprev];
    dsigma[*k] = d[jprev];
    idxp[*k] = jprev;

  L120:

//Count up the total number of the various types of columns, then
//form a permutation which positions the four column types into
//four groups of uniform structure (although one or more of these
//groups may be empty). */
    for (j = 1; j < 4; j++) {
	ctot[j - 1] = 0;
    }
    for (j = 2; j <= n; j++) {
	ct = coltyp[j];
	++ctot[ct - 1];
    }

//PSM(*) = Position in SubMatrix (of types 1 through 4)
    psm[0] = 2;
    psm[1] = ctot[0] + 2;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];

//Fill out the IDXC array so that the permutation which it induces
//will place all type-1 columns first, all type-2 columns next,
//then all type-3's, and finally all type-4's, starting from the
//second column. This applies similarly to the rows of VT.
    for (j = 2; j <= n; j++) {
	jp = idxp[j];
	ct = coltyp[jp];
	idxc[psm[ct - 1]] = j;
	++psm[ct - 1];

    }
//Sort the singular values and corresponding singular vectors into
//DSIGMA, U2, and VT2 respectively.  The singular values/vectors
//which were not deflated go into the first K slots of DSIGMA, U2,
//and VT2 respectively, while those which were deflated go into the
//last N - K slots, except that the first column/row will be treated
//separately. */
    for (j = 2; j <= n; j++) {
	jp = idxp[j];
	dsigma[j] = d[jp];
	idxj = idxq[idx[idxp[idxc[j]]] + 1];
	if (idxj <= nlp1) {
	    idxj--;
	}
	Rcopy(n, &u[idxj * ldu + 1], 1, &u2[j * ldu2 + 1], 1);
	Rcopy(m, &vt[idxj + ldvt], ldvt, &vt2[j + ldvt2], ldvt2);

    }
//Determine DSIGMA(1), DSIGMA(2) and Z(1)
    dsigma[1] = Zero;
    hlftol = tol / Two;
    if (abs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z[1] = Rlapy2(z1, z[m]);
	if (z[1] <= tol) {
	    c = One;
	    s = Zero;
	    z[1] = tol;
	} else {
	    c = z1 / z[1];
	    s = z[m] / z[1];
	}
    } else {
	if (abs(z1) <= tol) {
	    z[1] = tol;
	} else {
	    z[1] = z1;
	}
    }

//Move the rest of the updating row to Z.
    Rcopy(*k - 1, &u2[ldu2 + 2], 1, &z[2], 1);

//Determine the first column of U2, the first row of VT2 and the
//last row of VT.
    Rlaset("A", n, 1, Zero, Zero, &u2[0], ldu2);
    u2[nlp1 + ldu2] = One;
    if (m > n) {
	for (i = 0; i < nlp1; i++) {
	    vt[m + i * ldvt] = -s * vt[nlp1 + i * ldvt];
	    vt2[i * ldvt2 + 1] = c * vt[nlp1 + i * ldvt];

	}
	for (i = nlp2; i <= m; i++) {
	    vt2[i * ldvt2 + 1] = s * vt[m + i * ldvt];
	    vt[m + i * ldvt] = c * vt[m + i * ldvt];

	}
    } else {
	Rcopy(m, &vt[nlp1 + ldvt], ldvt, &vt2[ldvt2 + 1], ldvt2);
    }
    if (m > n) {
	Rcopy(m, &vt[m + ldvt], ldvt, &vt2[m + ldvt2], ldvt2);
    }
//The deflated singular values and their corresponding vectors go
//into the back of D, U, and V respectively.
    if (n > *k) {
	Rcopy(n - (*k), &dsigma[(*k) + 1], 1, &d[(*k) + 1], 1);
	Rlacpy("A", n, n - (*k), &u2[((*k) + 1) * ldu2 + 1], ldu2, &u[((*k) + 1) * ldu + 1], ldu);
	Rlacpy("A", n - (*k), m, &vt2[(*k) + 1 + ldvt2], ldvt2, &vt[(*k) + 1 + ldvt], ldvt);
    }
//Copy CTOT into COLTYP for referencing in DLASD3.
    for (j = 0; j < 4; j++) {
	coltyp[j] = ctot[j - 1];
    }
    return;
}
