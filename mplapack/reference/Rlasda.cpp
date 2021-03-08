/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasda.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasda(INTEGER icompq, INTEGER smlsiz, INTEGER n, INTEGER sqre, REAL * d, REAL * e,
       REAL * u, INTEGER ldu, REAL * vt, INTEGER * k, REAL * difl,
       REAL * difr, REAL * z, REAL * poles, INTEGER * givptr,
       INTEGER * givcol, INTEGER ldgcol, INTEGER * perm, REAL * givnum, REAL * c, REAL * s, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, m, i1, ic, lf, nd = 0, ll, nl, vf, nr, vl, im1, ncc, nlf, nrf, vfi, iwk, lvl, nru, ndb1, nlp1, lvl2, nrp1;
    REAL beta;
    INTEGER idxq, nlvl = 0;
    REAL alpha;
    INTEGER inode, ndiml, ndimr, idxqi, itemp;
    INTEGER nwork1, nwork2, smlszp, sqrei, vli;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (smlsiz < 3) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (sqre < 0 || sqre > 1) {
	*info = -4;
    } else if (ldu < n + sqre) {
	*info = -8;
    } else if (ldgcol < n) {
	*info = -17;
    }
    if (*info != 0) {
	Mxerbla("Rlasda", -(*info));
	return;
    }
    m = n + sqre;
//If the input matrix is too small, call DLASDQ to find the SVD.
    if (n <= smlsiz) {
	if (icompq == 0) {
	    Rlasdq("U", sqre, n, 0, 0, 0, &d[0], &e[0], &vt[0], ldu, &u[0], ldu, &u[0], ldu, &work[0], info);
	} else {
	    Rlasdq("U", sqre, n, m, n, 0, &d[0], &e[0], &vt[0], ldu, &u[0], ldu, &u[0], ldu, &work[0], info);
	}
	return;
    }
//Book-keeping and  set up the computation tree.
    inode = 1;
    ndiml = inode + n;
    ndimr = ndiml + n;
    idxq = ndimr + n;
    iwk = idxq + n;

    ncc = 0;
    nru = 0;

    smlszp = smlsiz + 1;
    vf = 1;
    vl = vf + m;
    nwork1 = vl + m;
    nwork2 = nwork1 + smlszp * smlszp;

    Rlasdt(n, nlvl, nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);

//for the nodes on bottom level of the tree, solve
//their subproblems by DLASDQ.
    ndb1 = (nd + 1) / 2;
    for (i = ndb1; i <= nd; i++) {
//IC : center row of each node
//NL : number of rows of left  subproblem
//NR : number of rows of right subproblem
//NLF: starting row of the left   subproblem
//NRF: starting row of the right  subproblem
	i1 = i - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i1];
	nlf = ic - nl;
	nrf = ic + 1;
	idxqi = idxq + nlf - 2;
	vfi = vf + nlf - 1;
	vli = vl + nlf - 1;
	sqrei = 0;

	if (icompq == 0) {
	    Rlaset("A", nlp1, nlp1, Zero, One, &work[nwork1], smlszp);
	    Rlasdq("U", sqrei, nl, nlp1, nru, ncc, &d[nlf], &e[nlf], &work[nwork1], smlszp, &work[nwork2], nl, &work[nwork2], nl, &work[nwork2], info);
	    itemp = nwork1 + nl * smlszp;
	    Rcopy(nlp1, &work[nwork1], 1, &work[vfi], 1);
	    Rcopy(nlp1, &work[itemp], 1, &work[vli], 1);
	} else {
	    Rlaset("A", nl, nl, Zero, One, &u[nlf + ldu], ldu);
	    Rlaset("A", nlp1, nlp1, Zero, One, &vt[nlf + ldu], ldu);
	    Rlasdq("U", sqrei, nl, nlp1, nl, ncc, &d[nlf], &e[nlf], &vt[nlf + ldu], ldu, &u[nlf + ldu], ldu, &u[nlf + ldu], ldu, &work[nwork1], info);
	    Rcopy(nlp1, &vt[nlf + ldu], 1, &work[vfi], 1);
	    Rcopy(nlp1, &vt[nlf + nlp1 * ldu], 1, &work[vli], 1);
	}
	if (*info != 0) {
	    return;
	}
	for (j = 0; j < nl; j++) {
	    iwork[idxqi + j] = j;
	}
	if (i == nd && sqre == 0) {
	    sqrei = 0;
	} else {
	    sqrei = 0;
	}
	idxqi += nlp1;
	vfi += nlp1;
	vli += nlp1;
	nrp1 = nr + sqrei;
	if (icompq == 0) {
	    Rlaset("A", nrp1, nrp1, Zero, One, &work[nwork1], smlszp);
	    Rlasdq("U", sqrei, nr, nrp1, nru, ncc, &d[nrf], &e[nrf], &work[nwork1], smlszp, &work[nwork2], nr, &work[nwork2], nr, &work[nwork2], info);
	    itemp = nwork1 + (nrp1 - 1) * smlszp;
	    Rcopy(nrp1, &work[nwork1], 1, &work[vfi], 1);
	    Rcopy(nrp1, &work[itemp], 1, &work[vli], 1);
	} else {
	    Rlaset("A", nr, nr, Zero, One, &u[nrf + ldu], ldu);
	    Rlaset("A", nrp1, nrp1, Zero, One, &vt[nrf + ldu], ldu);
	    Rlasdq("U", sqrei, nr, nrp1, nr, ncc, &d[nrf], &e[nrf], &vt[nrf + ldu], ldu, &u[nrf + ldu], ldu, &u[nrf + ldu], ldu, &work[nwork1], info);
	    Rcopy(nrp1, &vt[nrf + ldu], 1, &work[vfi], 1);
	    Rcopy(nrp1, &vt[nrf + nrp1 * ldu], 1, &work[vli], 1);
	}
	if (*info != 0) {
	    return;
	}
	for (j = 0; j < nr; j++) {
	    iwork[idxqi + j] = j;
	}
    }
//Now conquer each subproblem bottom-up.
    j = 2 ^ nlvl;
    for (lvl = nlvl; lvl >= 1; lvl--) {
	lvl2 = lvl * 2 - 1;
//Find the first node LF and last node LL on
//the current level LVL.
	if (lvl == 1) {
	    lf = 1;
	    ll = 0;
	} else {
	    lf = 2 ^ (lvl - 1);
	    ll = lf * 2 - 1;
	}
	for (i = lf; i <= ll; i++) {
	    im1 = i - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    nrf = ic + 1;
	    if (i == ll) {
		sqrei = sqre;
	    } else {
		sqrei = 0;
	    }
	    vfi = vf + nlf - 1;
	    vli = vl + nlf - 1;
	    idxqi = idxq + nlf - 1;
	    alpha = d[ic];
	    beta = e[ic];
	    if (icompq == 0) {
		Rlasd6(icompq, nl, nr, sqrei, &d[nlf], &work[vfi], &work[vli],
		       &alpha, &beta, &iwork[idxqi], &perm[0], &givptr[1],
		       &givcol[0], ldgcol, &givnum[0], ldu, &poles[0], &difl[0], &difr[0], &z[0], k[0], &c[1], &s[1], &work[nwork1], &iwork[iwk], info);
	    } else {
		j--;
		Rlasd6(icompq, nl, nr, sqrei, &d[nlf], &work[vfi], &work[vli],
		       &alpha, &beta, &iwork[idxqi], &perm[nlf + lvl * ldgcol],
		       &givptr[j], &givcol[nlf + lvl2 * ldgcol], ldgcol,
		       &givnum[nlf + lvl2 * ldgcol], ldu,
		       &poles[nlf + lvl2 * ldu], &difl[nlf + lvl * ldu], &difr[nlf + lvl2 * ldu], &z[nlf + lvl * ldu], k[j], &c[j], &s[j], &work[nwork1], &iwork[iwk], info);
	    }
	    if (*info != 0) {
		return;
	    }

	}
    }
    return;
}
