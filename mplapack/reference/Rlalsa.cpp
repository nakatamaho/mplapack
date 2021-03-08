/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlalsa.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlalsa(INTEGER icompq, INTEGER smlsiz, INTEGER n, INTEGER nrhs, REAL * B, INTEGER ldb,
       REAL * bx, INTEGER ldbx, REAL * u, INTEGER ldu, REAL * vt, INTEGER * k,
       REAL * difl, REAL * difr, REAL * z, REAL * poles,
       INTEGER * givptr, INTEGER * givcol, INTEGER ldgcol, INTEGER * perm, REAL * givnum, REAL * c, REAL * s, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, i1 = 0, ic, lf, nd = 0, ll, nl, nr, im1, nlf, nrf, lvl, ndb1, nlp1, lvl2, nrp1, nlvl = 0, sqre;
    INTEGER inode, ndiml, ndimr;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (smlsiz < 3) {
	*info = -2;
    } else if (n < smlsiz) {
	*info = -3;
    } else if (nrhs < 1) {
	*info = -4;
    } else if (ldb < n) {
	*info = -6;
    } else if (ldbx < n) {
	*info = -8;
    } else if (ldu < n) {
	*info = -10;
    } else if (ldgcol < n) {
	*info = -19;
    }
    if (*info != 0) {
	Mxerbla("Rlalsa", -(*info));
	return;
    }
//Book-keeping and  setting up the computation tree.
    inode = 1;
    ndiml = inode + n;
    ndimr = ndiml + n;
    Rlasdt(n, nlvl, nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
//The following code applies back the left singular vector factors.
//For applying back the right singular vector factors, go to 50
    if (icompq == 1) {
	goto L50;
    }
//The nodes on the bottom level of the tree were solved
//by DLASDQ. The corresponding left and right singular vector
//matrices are in explicit form. First apply back the left
//singular vector matrices.
    ndb1 = (nd + 1) / 2;
    for (i = ndb1; i <= nd; i++) {
//IC : center row of each node
//NL : number of rows of left  subproblem
//NR : number of rows of right subproblem
//NLF: starting row of the left   subproblem
//NRF: starting row of the right  subproblem
	ic = iwork[inode + i - 1];
	nl = iwork[ndiml + i - 1];
	nr = iwork[ndimr + i - 1];
	nlf = ic - nl;
	nrf = ic + 1;
	Rgemm("T", "N", nl, nrhs, nl, One, &u[nlf + ldu], ldu, &B[nlf + ldb], ldb, Zero, &bx[nlf + ldbx], ldbx);
	Rgemm("T", "N", nr, nrhs, nr, One, &u[nrf + ldu], ldu, &B[nrf + ldb], ldb, Zero, &bx[nrf + ldbx], ldbx);
    }
//Next copy the rows of B that correspond to unchanged rows
//in the bidiagonal matrix to BX.
    for (i = 0; i < nd; i++) {
	ic = iwork[inode + i - 1];
	Rcopy(nrhs, &B[ic + ldb], ldb, &bx[ic + ldbx], ldbx);
    }
//Finally go through the left singular vector matrices of all
//the other subproblems bottom-up on the tree.
    j = 2 ^ nlvl;
    sqre = 0;
    for (lvl = nlvl; lvl >= 1; lvl--) {
	lvl2 = (lvl * 2) - 1;
/*        find the first node LF and last node LL on */
/*        the current level LVL */

	if (lvl == 1) {
	    lf = 1;
	    ll = 0;
	} else {
	    lf = 2 ^ (lvl - 1);
	    ll = (lf * 2) - 1;
	}
	for (i = lf; i <= i1; i++) {
	    im1 = i - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    nrf = ic + 1;
	    j--;
	    //      Rlals0(icompq, nl, nr, sqre, nrhs, &bx[nlf + ldbx], ldbx, &B[nlf + ldb], ldb, &perm[nlf + lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * givcol_dim1],ldgcol, &givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * polelds], &difl[nlf + lvl * lddifl], &difr[nlf + lvl2 * lddifr],&z[nlf + lvl * ldz], &k[j], &c[j], &s[j], &work[0], info);
	}
    }
    goto L90;
//ICOMPQ = 1: applying back the right singular vector factors.
  L50:
//First now go through the right singular vector matrices of all
//the tree nodes top-down.
    j = 0;
    for (lvl = 0; lvl <= nlvl; lvl++) {
	lvl2 = (lvl * 2) - 1;
//Find the first node LF and last node LL on
//the current level LVL.
	if (lvl == 1) {
	    lf = 1;
	    ll = 0;
	} else {
	    lf = 2 ^ (lvl - 1);
	    ll = (lf * 2) - 1;
	}
	for (i = ll; i >= lf; i--) {
	    im1 = i - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    nrf = ic + 1;
	    if (i == ll) {
		sqre = 0;
	    } else {
		sqre = 1;
	    }
	    j++;
	    //Rlals0(icompq, &nl, &nr, &sqre, nrhs, &B[nlf + ldb], ldb, &bx[nlf + ldbx], ldbx, &perm[nlf + lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * givcol_dim1],ldgcol, &givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * polelds], &difl[nlf + lvl * lddifl], &difr[nlf + lvl2 * lddifr], &z[nlf + lvl * ldz], &k[j], &c[j], &s[j], &work[0], info);
	}
    }
//The nodes on the bottom level of the tree were solved
//by DLASDQ. The corresponding right singular vector
//matrices are in explicit form. Apply them back.
    ndb1 = (nd + 1) / 2;
    for (i = ndb1; i <= nd; i++) {
	i1 = i - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nr = iwork[ndimr + i1];
	nlp1 = nl + 1;
	if (i == nd) {
	    nrp1 = nr;
	} else {
	    nrp1 = nr + 1;
	}
	nlf = ic - nl;
	nrf = ic + 1;
	Rgemm("T", "N", nlp1, nrhs, nlp1, One, &vt[nlf + ldu], ldu, &B[nlf + ldb], ldb, Zero, &bx[nlf + ldbx], ldbx);
	Rgemm("T", "N", nrp1, nrhs, nrp1, One, &vt[nrf + ldu], ldu, &B[nrf + ldb], ldb, Zero, &bx[nrf + ldbx], ldbx);

    }
  L90:
    return;
}
