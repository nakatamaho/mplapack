/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd7.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasd7(INTEGER icompq, INTEGER nl, INTEGER nr, INTEGER sqre, INTEGER k, REAL * d,
       REAL * z, REAL * zw, REAL * vf, REAL * vfw,
       REAL * vl, REAL * vlw, REAL alpha, REAL beta,
       REAL * dsigma, INTEGER * idx, INTEGER * idxp, INTEGER * idxq, INTEGER * perm, INTEGER * givptr,
       INTEGER * givcol, INTEGER ldgcol, REAL * givnum, INTEGER ldgnum, REAL * c, REAL * s, INTEGER * info)
{
    INTEGER i, j, m, n, k2;
    REAL z1;
    INTEGER jp;
    REAL eps, tau, tol;
    INTEGER nlp1, nlp2, idxi, idxj;
    INTEGER idxjp;
    INTEGER jprev = 0;
    REAL hlftol;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Eight = 8.0;
    REAL mtemp1;

    *info = 0;
    n = nl + nr + 1;
    m = n + sqre;

    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (nl < 1) {
	*info = -2;
    } else if (nr < 1) {
	*info = -3;
    } else if (sqre < 0 || sqre > 1) {
	*info = -4;
    } else if (ldgcol < n) {
	*info = -22;
    } else if (ldgnum < n) {
	*info = -24;
    }
    if (*info != 0) {
	Mxerbla("Rlasd7", -(*info));
	return;
    }

    nlp1 = nl + 1;
    nlp2 = nl + 2;
    if (icompq == 1) {
	*givptr = 0;
    }
//Generate the first part of the vector Z and move the singular
//values in the first part of D one position backward.
    z1 = alpha * vl[nlp1];
    vl[nlp1] = Zero;
    tau = vf[nlp1];
    for (i = nl; i >= 1; i--) {
	z[i + 1] = alpha * vl[i];
	vl[i] = Zero;
	vf[i + 1] = vf[i];
	d[i + 1] = d[i];
	idxq[i + 1] = idxq[i] + 1;
    }
    vf[1] = tau;
//Generate the second part of the vector Z.
    for (i = nlp2; i <= m; i++) {
	z[i] = beta * vf[i];
	vf[i] = Zero;
    }
//Sort the singular values into increasing order
    for (i = nlp2; i <= n; i++) {
	idxq[i] = idxq[i] + nlp1;
    }
//DSIGMA, IDXC, IDXC, and ZW are used as storage space.
    for (i = 1; i < n; i++) {
	dsigma[i] = d[idxq[i]];
	zw[i] = z[idxq[i]];
	vfw[i] = vf[idxq[i]];
	vlw[i] = vl[idxq[i]];
    }
    Rlamrg(nl, nr, &dsigma[2], 1, 1, &idx[2]);
    for (i = 1; i < n; i++) {
	idxi = idx[i] + 1;
	d[i] = dsigma[idxi];
	z[i] = zw[idxi];
	vf[i] = vfw[idxi];
	vl[i] = vlw[idxi];
    }
//Calculate the allowable deflation tolerence
    eps = Rlamch("Epsilon");
    tol = max(abs(alpha), abs(beta));
    mtemp1 = abs(d[n]);
    tol = eps * Eight * Eight * max(mtemp1, tol);
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
	    if (j == n) {
		goto L100;
	    }
	} else {
	    jprev = j;
	    goto L70;
	}
    }
  L70:
    j = jprev;
  L80:
    j++;
    if (j > n) {
	goto L90;
    }
    if (abs(z[j]) <= tol) {
//Deflate due to small z component.
	k2--;
	idxp[k2] = j;
    } else {
//Check if singular values are close enough to allow deflation.
	if (abs(d[j] - d[jprev]) <= tol) {
//Deflation is possible.
	    *s = z[jprev];
	    *c = z[j];
//Find sqrt(a**2+b**2) without overflow or
//destructive underflow.
	    tau = Rlapy2(*c, *s);
	    z[j] = tau;
	    z[jprev] = Zero;
	    *c = *c / tau;
	    *s = -(*s) / tau;
//Record the appropriate Givens rotation
	    if (icompq == 1) {
		++(*givptr);
		idxjp = idxq[idx[jprev] + 1];
		idxj = idxq[idx[j] + 1];
		if (idxjp <= nlp1) {
		    idxjp--;
		}
		if (idxj <= nlp1) {
		    idxj--;
		}
		givcol[*givptr + (ldgcol * 2)] = idxjp;
		givcol[*givptr + ldgcol] = idxj;
		givnum[*givptr + (ldgnum * 2)] = *c;
		givnum[*givptr + ldgnum] = *s;
	    }
	    Rrot(1, &vf[jprev], 1, &vf[j], 1, *c, *s);
	    Rrot(1, &vl[jprev], 1, &vl[j], 1, *c, *s);
	    k2--;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(k);
	    zw[k] = z[jprev];
	    dsigma[k] = d[jprev];
	    idxp[k] = jprev;
	    jprev = j;
	}
    }
    goto L80;
  L90:
//Record the last singular value.
    ++(k);
    zw[k] = z[jprev];
    dsigma[k] = d[jprev];
    idxp[k] = jprev;

  L100:

//Sort the singular values into DSIGMA. The singular values which
//were not deflated go into the first K slots of DSIGMA, except
//that DSIGMA(1) is treated separately.
    for (j = 2; j <= n; j++) {
	jp = idxp[j];
	dsigma[j] = d[jp];
	vfw[j] = vf[jp];
	vlw[j] = vl[jp];
    }
    if (icompq == 1) {
	for (j = 2; j <= n; j++) {
	    jp = idxp[j];
	    perm[j] = idxq[idx[jp] + 1];
	    if (perm[j] <= nlp1) {
		--perm[j];
	    }
	}
    }
//The deflated singular values go back into the last N - K slots of
//D.
    Rcopy(n - k, &dsigma[k + 1], 1, &d[k + 1], 1);
//Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
//VL(M).
    dsigma[1] = Zero;
    hlftol = tol / Two;
    if (abs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z[1] = Rlapy2(z1, z[m]);
	if (z[1] <= tol) {
	    *c = One;
	    *s = Zero;
	    z[1] = tol;
	} else {
	    *c = z1 / z[1];
	    *s = -z[m] / z[1];
	}
	Rrot(1, &vf[m], 1, &vf[1], 1, *c, *s);
	Rrot(1, &vl[m], 1, &vl[1], 1, *c, *s);
    } else {
	if (abs(z1) <= tol) {
	    z[1] = tol;
	} else {
	    z[1] = z1;
	}
    }
//Restore Z, VF, and VL.
    Rcopy(k - 1, &zw[2], 1, &z[2], 1);
    Rcopy(n - 1, &vfw[2], 1, &vf[2], 1);
    Rcopy(n - 1, &vlw[2], 1, &vl[2], 1);
    return;
}
