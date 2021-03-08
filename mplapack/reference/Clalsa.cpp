/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clalsa.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Clalsa(INTEGER icompq, INTEGER smlsiz, INTEGER n, INTEGER nrhs, COMPLEX * B, INTEGER ldb,
       COMPLEX * bx, INTEGER ldbx, REAL * u, INTEGER ldu, REAL * vt, INTEGER * k,
       REAL * difl, REAL * difr, REAL * z, REAL * poles,
       INTEGER * givptr, INTEGER * givcol, INTEGER ldgcol, INTEGER * perm, REAL * givnum, REAL * c, REAL * s, REAL * rwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, i1, ic, lf, nd = 0, ll, nl, nr, im1, nlf, nrf, lvl, ndb1, nlp1, lvl2, nrp1, jcol, nlvl = 0, sqre, jrow, jimag;
    INTEGER jreal, inode, ndiml, ndimr;
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
	Mxerbla("Clalsa", -(*info));
	return;
    }
//Book-keeping and  setting up the computation tree.
    inode = 1;
    ndiml = inode + n;
    ndimr = ndiml + n;
    Rlasdt(n, nlvl, nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);

//The following code applies back the left singular vector factors.
//For applying back the right singular vector factors, go to 170
    if (icompq == 1) {
	goto L170;
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
	i1 = i - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nr = iwork[ndimr + i1];
	nlf = ic - nl;
	nrf = ic + 1;

// Since B and BX are complex, the following call to DGEMM
// is performed in two steps (real and imaginary parts).
// CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
//              B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
	j = nl * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nlf; jrow <= nlf + nl - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].real();
	    }

	}
	Rgemm("T", "N", nl, nrhs, nl, One, &u[nlf + ldu], ldu, &rwork[(nl * nrhs * 2) + 1], nl, Zero, &rwork[1], nl);
	j = nl * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nlf; jrow <= nlf + nl - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].imag();
	    }
	}
	Rgemm("T", "N", nl, nrhs, nl, One, &u[nlf + ldu], ldu, &rwork[(nl * nrhs * 2) + 1], nl, Zero, &rwork[nl * nrhs + 1], nl);
	jreal = 0;
	jimag = nl * nrhs;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nlf; jrow <= nlf + nl - 1; jrow++) {
		jreal++;
		jimag++;
		bx[jrow + jcol * ldbx] = rwork[jreal];
	    }
	}
//Since B and BX are complex, the following call to DGEMM
//is performed in two steps (real and imaginary parts).
//CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
//            B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
	j = nr * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nrf; jrow <= nrf + nr - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].real();
	    }
	}
	Rgemm("T", "N", nr, nrhs, nr, One, &u[nrf + ldu], ldu, &rwork[(nr * nrhs * 2) + 1], nr, Zero, &rwork[1], nr);
	j = nr * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nrf; jrow <= nrf + nr - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].imag();
	    }
	}
	Rgemm("T", "N", nr, nrhs, nr, One, &u[nrf + ldu], ldu, &rwork[(nr * nrhs * 2) + 1], nr, Zero, &rwork[nr * nrhs + 1], nr);
	jreal = 0;
	jimag = nr * nrhs;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nrf; jrow <= nrf + nr - 1; jrow++) {
		jreal++;
		jimag++;
		bx[jrow + jcol * ldbx] = rwork[jreal];
	    }
	}
    }
//Next copy the rows of B that correspond to unchanged rows
//in the bidiagonal matrix to BX.
    for (i = 0; i < nd; i++) {
	ic = iwork[inode + i - 1];
	Ccopy(nrhs, &B[ic + ldb], ldb, &bx[ic + ldbx], ldbx);
    }
//Finally go through the left singular vector matrices of all
//the other subproblems bottom-up on the tree.
    j = 2 ^ nlvl;
    sqre = 0;
    for (lvl = nlvl; lvl >= 1; lvl--) {
	lvl2 = (lvl * 2) - 1;
//find the first node LF and last node LL on
//the current level LVL
	if (lvl == 1) {
	    lf = 1;
	    ll = 0;
	} else {
	    lf = 2 ^ (lvl - 1);
	    ll = (lf * 2) - 1;
	}
	for (i = lf; i <= ll; i++) {
	    im1 = i - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    nrf = ic + 1;
	    j--;
	    Clals0(icompq, nl, nr, sqre, nrhs, &bx[nlf + ldbx], ldbx,
		   &B[nlf + ldb], ldb, &perm[nlf + lvl * ldgcol], givptr[j],
		   &givcol[nlf + lvl2 * ldgcol], ldgcol,
		   &givnum[nlf + lvl2 * ldu], ldu,
		   &poles[nlf + lvl2 * ldu], &difl[nlf + lvl * ldu], &difr[nlf + lvl2 * ldu], &z[nlf + lvl * ldu], k[j], c[j], s[j], &rwork[1], info);
	}
    }
    goto L330;
//ICOMPQ = 1: applying back the right singular vector factors.
  L170:
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
	    Clals0(icompq, nl, nr, sqre, nrhs, &B[nlf + ldb], ldb,
		   &bx[nlf + ldbx], ldbx, &perm[nlf + lvl * ldgcol],
		   givptr[j], &givcol[nlf + lvl2 * ldgcol], ldgcol,
		   &givnum[nlf + lvl2 * ldu], ldu,
		   &poles[nlf + lvl2 * ldu], &difl[nlf + lvl * ldu], &difr[nlf + lvl2 * ldu], &z[nlf + lvl * ldu], k[j], c[j], s[j], &rwork[1], info);
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
//Since B and BX are complex, the following call to DGEMM is
//performed in two steps (real and imaginary parts).
//CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
//            B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
	j = nlp1 * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nlf; jrow <= nlf + nlp1 - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].real();
	    }
	}
	Rgemm("T", "N", nlp1, nrhs, nlp1, One, &vt[nlf + ldu], ldu, &rwork[(nlp1 * nrhs * 2) + 1], nlp1, Zero, &rwork[1], nlp1);
	j = nlp1 * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nlf; jrow <= nlf + nlp1 - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].imag();
	    }
	}
	Rgemm("T", "N", nlp1, nrhs, nlp1, One, &vt[nlf + ldu], ldu, &rwork[(nlp1 * nrhs * 2) + 1], nlp1, Zero, &rwork[nlp1 * nrhs + 1], nlp1);
	jreal = 0;
	jimag = nlp1 * nrhs;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nlf; jrow <= nlf + nlp1 - 1; jrow++) {
		jreal++;
		jimag++;
		bx[jrow + jcol * ldbx] = (COMPLEX) rwork[jreal];
	    }
	}
//Since B and BX are complex, the following call to DGEMM is
//performed in two steps (real and imaginary parts).
//CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
//            B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
	j = nrp1 * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nrf; jrow <= nrf + nrp1 - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].real();
	    }
	}
	Rgemm("T", "N", nrp1, nrhs, nrp1, One, &vt[nrf + ldu], ldu, &rwork[(nrp1 * nrhs * 2) + 1], nrp1, Zero, &rwork[1], nrp1);
	j = nrp1 * nrhs * 2;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nrf; jrow <= nrf + nrp1 - 1; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].imag();
	    }
	}
	Rgemm("T", "N", nrp1, nrhs, nrp1, One, &vt[nrf + ldu], ldu, &rwork[(nrp1 * nrhs * 2) + 1], nrp1, Zero, &rwork[nrp1 * nrhs + 1], nrp1);
	jreal = 0;
	jimag = nrp1 * nrhs;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = nrf; jrow <= nrf + nrp1 - 1; jrow++) {
		jreal++;
		jimag++;
		bx[jrow + jcol * ldbx] = rwork[jreal];
	    }
	}
    }
  L330:
    return;
}
