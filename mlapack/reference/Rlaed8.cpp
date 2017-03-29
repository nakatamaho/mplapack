/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed8.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaed8(INTEGER icompq, INTEGER * k, INTEGER n, INTEGER qsiz, REAL * d, REAL * q,
       INTEGER ldq, INTEGER * indxq, REAL * rho, INTEGER cutpnt, REAL * z,
       REAL * dlamda, REAL * q2, INTEGER ldq2, REAL * w, INTEGER * perm, INTEGER * givptr, INTEGER * givcol, REAL * givnum, INTEGER * indxp, INTEGER * indx, INTEGER * info)
{

    REAL c;
    INTEGER i, j;
    REAL s, t;
    INTEGER k2, n1, n2, jp, n1p1;
    REAL eps, tau, tol;
    INTEGER jlam = 0, imax, jmax;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Eight = 8.0;

    *info = 0;
    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (n < 0) {
	*info = -3;
    } else if (icompq == 1 && qsiz < n) {
	*info = -4;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (cutpnt < max((INTEGER) 1, n) || cutpnt > n) {
	*info = -10;
    } else if (ldq2 < max((INTEGER) 1, n)) {
	*info = -14;
    }
    if (*info != 0) {
	Mxerbla("Rlaed8", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    n1 = cutpnt;
    n2 = n - n1;
    n1p1 = n1 + 1;
    if (*rho < Zero) {
	Rscal(n2, -One, &z[n1p1], 1);
    }
//Normalize z so that norm(z) = 1
    t = One / sqrt(Two);
    for (j = 0; j < n; j++) {
	indx[j] = j;
    }
    Rscal(n, t, &z[0], 1);
    *rho = abs(*rho * Two);
//Sort the eigenvalues into increasing order
    for (i = cutpnt + 1; i <= n; i++) {
	indxq[i] += cutpnt;
    }
    for (i = 0; i < n; i++) {
	dlamda[i] = d[indxq[i]];
	w[i] = z[indxq[i]];

    }
    i = 1;
    j = cutpnt + 1;
    Rlamrg(n1, n2, &dlamda[1], 1, 1, &indx[0]);
    for (i = 0; i < n; i++) {
	d[i] = dlamda[indx[i]];
	z[i] = w[indx[i]];
    }
//Calculate the allowable deflation tolerence
    imax = iRamax(n, &z[1], 1);
    jmax = iRamax(n, &d[0], 1);
    eps = Rlamch("Epsilon");
    tol = eps * Eight * abs(d[jmax]);
//If the rank-1 modifier is small enough, no more needs to be done
//except to reorganize Q so that its columns correspond with the
//elements in D.
    if (*rho * abs(z[imax]) <= tol) {
	k = 0;
	if (icompq == 0) {
	    for (j = 0; j < n; j++) {
		perm[j] = indxq[indx[j]];

	    }
	} else {
	    for (j = 0; j < n; j++) {
		perm[j] = indxq[indx[j]];
		Rcopy(qsiz, &q[perm[j] * ldq + 1], 1, &q2[j * ldq2 + 1], 1);

	    }
	    Rlacpy("A", qsiz, n, &q2[ldq2 + 1], ldq2, &q[ldq + 1], ldq);
	}
	return;
    }
//If there are multiple eigenvalues then the problem deflates.  Here
//the number of equal eigenvalues are found.  As each equal
//eigenvalue is found, an elementary reflector is computed to rotate
//the corresponding eigensubspace so that the corresponding
//components of Z are zero in this new basis.

    k = 0;
    *givptr = 0;
    k2 = n + 1;
    for (j = 0; j < n; j++) {
	if (*rho * abs(z[j]) <= tol) {
//Deflate due to small z component.
	    k2--;
	    indxp[k2] = j;
	    if (j == n) {
		goto L110;
	    }
	} else {
	    jlam = j;
	    goto L80;
	}
    }
  L80:
    j++;
    if (j > n) {
	goto L100;
    }
    if (*rho * abs(z[j]) <= tol) {
//Deflate due to small z component.
	k2--;
	indxp[k2] = j;
    } else {
//Check if eigenvalues are close enough to allow deflation.
	s = z[jlam];
	c = z[j];
//Find sqrt(a**2+b**2) without overflow or
//destructive underflow.
	tau = Rlapy2(c, s);
	t = d[j] - d[jlam];
	c /= tau;
	s = -s / tau;
	if (abs(t * c * s) <= tol) {
//Deflation is possible.
	    z[j] = tau;
	    z[jlam] = Zero;
//Record the appropriate Givens rotation
	    ++(*givptr);
	    givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
	    givcol[(*givptr << 1) + 2] = indxq[indx[j]];
	    givnum[(*givptr << 1) + 1] = c;
	    givnum[(*givptr << 1) + 2] = s;
	    if (icompq == 1) {
		Rrot(qsiz, &q[indxq[indx[jlam]] * ldq + 1], 1, &q[indxq[indx[j]] * ldq + 1], 1, c, s);
	    }
	    t = d[jlam] * c * c + d[j] * s * s;
	    d[j] = d[jlam] * s * s + d[j] * c * c;
	    d[jlam] = t;
	    k2--;
	    i = 1;
	  L90:
	    if (k2 + i <= n) {
		if (d[jlam] < d[indxp[k2 + i]]) {
		    indxp[k2 + i - 1] = indxp[k2 + i];
		    indxp[k2 + i] = jlam;
		    i++;
		    goto L90;
		} else {
		    indxp[k2 + i - 1] = jlam;
		}
	    } else {
		indxp[k2 + i - 1] = jlam;
	    }
	    jlam = j;
	} else {
	    ++(*k);
	    w[*k] = z[jlam];
	    dlamda[*k] = d[jlam];
	    indxp[*k] = jlam;
	    jlam = j;
	}
    }
    goto L80;
  L100:

//Record the last eigenvalue.
    ++(*k);
    w[*k] = z[jlam];
    dlamda[*k] = d[jlam];
    indxp[*k] = jlam;

  L110:
//Sort the eigenvalues and corresponding eigenvectors into DLAMDA
//and Q2 respectively.  The eigenvalues/vectors which were not
//deflated go into the first K slots of DLAMDA and Q2 respectively,
//while those which were deflated go into the last N - K slots.
    if (icompq == 0) {
	for (j = 0; j < n; j++) {
	    jp = indxp[j];
	    dlamda[j] = d[jp];
	    perm[j] = indxq[indx[jp]];
	}
    } else {
	for (j = 0; j < n; j++) {
	    jp = indxp[j];
	    dlamda[j] = d[jp];
	    perm[j] = indxq[indx[jp]];
	    Rcopy(qsiz, &q[perm[j] * ldq + 1], 1, &q2[j * ldq2 + 1], 1);

	}
    }
//The deflated eigenvalues and their corresponding vectors go back
//into the last N - K slots of D and Q respectively.
    if (*k < n) {
	if (icompq == 0) {
	    Rcopy(n - *k, &dlamda[*k + 1], 1, &d[*k + 1], 1);
	} else {
	    Rcopy(n - *k, &dlamda[*k + 1], 1, &d[*k + 1], 1);
	    Rlacpy("A", qsiz, n - *k, &q2[(*k + 1) * ldq2 + 1], ldq2, &q[(*k + 1) * ldq + 1], ldq);
	}
    }
    return;
}
