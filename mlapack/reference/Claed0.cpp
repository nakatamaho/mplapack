/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claed0.cpp,v 1.13 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Claed0(INTEGER qsiz, INTEGER n, REAL * d, REAL * e, COMPLEX * q, INTEGER ldq, COMPLEX * qstore, INTEGER ldqs, REAL * rwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, k, ll, iq, lgn, msd2, smm1, spm1, spm2;
    REAL temp;
    INTEGER curr, iperm;
    INTEGER indxq, iwrem, iqptr, tlvls;
    INTEGER igivcl;
    INTEGER igivnm, submat, curprb = 0, subpbs, igivpt;
    INTEGER curlvl, matsiz, iprmpt, smlsiz;
    REAL mtemp1;

//Warning:      N could be as big as QSIZ!
//Test the input parameters.
    *info = 0;
    if (qsiz < max((INTEGER) 0, n)) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldqs < max((INTEGER) 1, n)) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Claed0", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }

    smlsiz = iMlaenv(9, "Claed0", " ", 0, 0, 0, 0);
//Determine the size and placement of the submatrices, and save in
//the leading elements of IWORK.
    iwork[1] = n;
    subpbs = 1;
    tlvls = 0;
  L10:
    if (iwork[subpbs] > smlsiz) {
	for (j = subpbs; j >= 1; j--) {
	    iwork[j * 2] = (iwork[j] + 1) / 2;
	    iwork[(j << 1) - 1] = iwork[j] / 2;

	}
	++tlvls;
	subpbs <<= 1;
	goto L10;
    }
    for (j = 2; j <= subpbs; j++) {
	iwork[j] += iwork[j - 1];
    }
//Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
//using rank-1 modifications (cuts).
    spm1 = subpbs - 1;
    for (i = 0; i < spm1; i++) {
	submat = iwork[i] + 1;
	smm1 = submat - 1;
	d[smm1] = d[smm1] - abs(e[smm1]);
	d[submat] = d[submat] - abs(e[smm1]);

    }
    indxq = (n << 2) + 3;
//Set up workspaces for eigenvalues only/accumulate new vectors
//routine
    temp = (double) n;
    mtemp1 = log(temp);
    lgn = (INTEGER) cast2double(mtemp1);
    if ((2 ^ lgn) < n) {
	lgn++;
    }
    if ((2 ^ lgn) < n) {
	lgn++;
    }
    iprmpt = indxq + n + 1;
    iperm = iprmpt + n * lgn;
    iqptr = iperm + n * lgn;
    igivpt = iqptr + n + 2;
    igivcl = igivpt + n * lgn;
    igivnm = 1;
    iq = igivnm + (n << 1) * lgn;

    iwrem = iq + n * n + 1;
//Initialize pointers
    for (i = 0; i < subpbs; i++) {
	iwork[iprmpt + i] = 1;
	iwork[igivpt + i] = 1;
    }
    iwork[iqptr] = 1;
//Solve each submatrix eigenproblem at the bottom of the divide and
//conquer tree.
    curr = 0;
    for (i = 0; i < spm1; i++) {
	if (i == 0) {
	    submat = 1;
	    matsiz = iwork[1];
	} else {
	    submat = iwork[i] + 1;
	    matsiz = iwork[i + 1] - iwork[i];
	}
	ll = iq - 1 + iwork[iqptr + curr];
	Rsteqr("I", matsiz, &d[submat], &e[submat], &rwork[ll], matsiz, &rwork[1], info);
	Clacrm(qsiz, matsiz, &q[submat * ldq + 1], ldq, &rwork[ll], matsiz, &qstore[submat * ldq + 1], ldqs, &rwork[iwrem]);
	iwork[iqptr + curr + 1] = iwork[iqptr + curr] + matsiz * matsiz;
	++curr;
	if (*info > 0) {
	    *info = submat * (n + 1) + submat + matsiz - 1;
	    return;
	}
	k = 0;
	for (j = submat; j <= iwork[i + 1]; j++) {
	    iwork[indxq + j] = k;
	    k++;
	}
    }
//Successively merge eigensystems of adjacent submatrices
//into eigensystem for the corresponding larger matrix.
    curlvl = 0;
    while (subpbs > 1) {
	spm2 = subpbs - 2;
	for (i = 0; i < spm2; i += 2) {
	    if (i == 0) {
		submat = 1;
		matsiz = iwork[2];
		msd2 = iwork[1];
		curprb = 0;
	    } else {
		submat = iwork[i] + 1;
		matsiz = iwork[i + 2] - iwork[i];
		msd2 = matsiz / 2;
		++curprb;
	    }

//Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
//into an eigensystem of size MATSIZ.  ZLAED7 handles the case
//when the eigenvectors of a full or band Hermitian matrix (which
//was reduced to tridiagonal form) are desired.
//I am free to use Q as a valuable working space until Loop 150
	    Claed7(matsiz, msd2, qsiz, tlvls, curlvl, curprb, &d[submat],
		   &qstore[submat * ldq + 1], ldqs, e[submat + msd2 - 1],
		   &iwork[indxq + submat], &rwork[iq], &iwork[iqptr],
		   &iwork[iprmpt], &iwork[iperm], &iwork[igivpt], &iwork[igivcl], &rwork[igivnm], &q[submat * ldq + 1], &rwork[iwrem], &iwork[subpbs + 1], info);
	    if (*info > 0) {
		*info = submat * (n + 1) + submat + matsiz - 1;
		return;
	    }
	    iwork[i / 2 + 1] = iwork[i + 2];

	}
	subpbs /= 2;
	++curlvl;
    }

//Re-merge the eigenvalues/vectors which were deflated at the final
//merge step.
    for (i = 0; i < n; i++) {
	j = iwork[indxq + i];
	rwork[i] = d[j];
	Ccopy(qsiz, &qstore[j * ldq + 1], 1, &q[i * ldq + 1], 1);

    }
    Rcopy(n, &rwork[1], 1, &d[0], 1);
    return;
}
