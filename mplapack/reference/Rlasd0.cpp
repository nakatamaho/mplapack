/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd0.cpp,v 1.10 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasd0(INTEGER n, INTEGER sqre, REAL * d, REAL * e, REAL * u, INTEGER ldu, REAL * vt, INTEGER ldvt, INTEGER smlsiz, INTEGER * iwork, REAL * work, INTEGER * info)
{
    INTEGER i, j, m, ic, lf, nd, ll, nl, nr, im1, ncc, nlf, nrf, iwk, lvl, ndb1, nlp1, nrp1;
    REAL beta;
    INTEGER idxq, nlvl = 0.0;
    REAL alpha;
    INTEGER inode, ndiml, idxqc, ndimr, itemp, sqrei;

    nd = 0.0;
//Test the input parameters.
    *info = 0;
    if (n < 0) {
	*info = -1;
    } else if (sqre < 0 || sqre > 1) {
	*info = -2;
    }
    m = n + sqre;
    if (ldu < n) {
	*info = -6;
    } else if (ldvt < m) {
	*info = -8;
    } else if (smlsiz < 3) {
	*info = -9;
    }
    if (*info != 0) {
	Mxerbla("Rlasd0", -(*info));
	return;
    }
//If the input matrix is too small, call DLASDQ to find the SVD.
    if (n <= smlsiz) {
	Rlasdq("U", sqre, n, m, n, 0, &d[0], &e[0], &vt[0], ldvt, &u[0], ldu, &u[0], ldu, &work[0], info);
	return;
    }
//Set up the computation tree.
    inode = 1;
    ndiml = inode + n;
    ndimr = ndiml + n;
    idxq = ndimr + n;
    iwk = idxq + n;
    Rlasdt(n, nlvl, nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
//For the nodes on bottom level of the tree, solve
//their subproblems by DLASDQ.
    ndb1 = (nd + 1) / 2;
    ncc = 0;
    for (i = ndb1; i <= nd; i++) {
//IC : center row of each node
//NL : number of rows of left  subproblem
//NR : number of rows of right subproblem
//NLF: starting row of the left   subproblem
//NRF: starting row of the right  subproblem
	ic = iwork[inode + i - 1];
	nl = iwork[ndiml + i - 1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i - 1];
	nrp1 = nr + 1;
	nlf = ic - nl;
	nrf = ic + 1;
	sqre = 0;

	Rlasdq("U", sqrei, nl, nlp1, nl, ncc, &d[nlf], &e[nlf], &vt[nlf + nlf * ldvt], ldvt, &u[nlf + nlf * ldu], ldu, &u[nlf + nlf * ldu], ldu, &work[0], info);
	if (*info != 0) {
	    return;
	}
	itemp = idxq + nlf - 2;
	for (j = 0; j < nl; j++) {
	    iwork[itemp + j] = j;
	}
	if (i == nd) {
	    sqrei = sqre;
	} else {
	    sqrei = 0;
	}
	nrp1 = nr + sqrei;
	Rlasdq("U", sqrei, nr, nrp1, nr, ncc, &d[nrf], &e[nrf], &vt[nrf + nrf * ldvt], ldvt, &u[nrf + nrf * ldu], ldu, &u[nrf + nrf * ldu], ldu, &work[0], info);
	if (*info != 0) {
	    return;
	}
	itemp = idxq + ic;
	for (j = 0; j < nr; j++) {
	    iwork[itemp + j - 1] = j;
	}
    }
//Now conquer each subproblem bottom-up.
    for (lvl = nlvl; lvl >= 1; lvl--) {
//Find the first node LF and last node LL on the
//current level LVL.
	if (lvl == 1) {
	    lf = 1;
	    ll = 0;
	} else {
	    lf = 2 ^ (lvl - 1);
	    ll = (lf << 1) - 1;
	}
	for (i = lf; i <= ll; i++) {
	    im1 = i - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    if (sqre == 0 && i == ll) {
		sqrei = sqre;
	    } else {
		sqrei = 0;
	    }
	    idxqc = idxq + nlf - 1;
	    alpha = d[ic];
	    beta = e[ic];
	    Rlasd1(nl, nr, &sqrei, &d[nlf], &alpha, &beta, &u[nlf + nlf * ldu], ldu, &vt[nlf + nlf * ldvt], ldvt, &iwork[idxqc], &iwork[iwk], &work[0], info);
	    if (*info != 0) {
		return;
	    }
	}
    }
    return;
}
