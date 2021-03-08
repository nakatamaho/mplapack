/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed0.cpp,v 1.14 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlaed0(INTEGER icompq, INTEGER qsiz, INTEGER n, REAL * d, REAL * e, REAL * q, INTEGER ldq, REAL * qstore, INTEGER ldqs, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, k, iq = 0, lgn, msd2, smm1, spm1, spm2;
    REAL temp;
    INTEGER curr;
    INTEGER iperm = 0;
    INTEGER indxq, iwrem = 0;
    INTEGER iqptr = 0;
    INTEGER tlvls;
    INTEGER igivcl = 0;
    INTEGER igivnm, submat, curprb = 0, subpbs, igivpt = 0;
    INTEGER curlvl, matsiz, iprmpt = 0, smlsiz;
    REAL mtemp1;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    if (icompq < 0 || icompq > 2) {
	*info = -1;
    } else if (icompq == 1 && qsiz < max((INTEGER) 0, n)) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldqs < max((INTEGER) 1, n)) {
	*info = -9;
    }
    if (*info != 0) {
	Mxerbla("Rlaed0", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    smlsiz = iMlaenv(9, "Rlaed", " ", 0, 0, 0, 0);
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
    if (icompq != 2) {
//Set up workspaces for eigenvalues only/accumulate new vectors
//routine
	mtemp1 = (double) n;
	temp = log(mtemp1);
	lgn = (INTEGER) cast2double(temp);
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
    }
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
	if (icompq == 2) {
	    Rsteqr("I", matsiz, &d[submat], &e[submat], &q[submat + submat * ldq], ldq, &work[0], info);
	    if (*info != 0) {
		goto L130;
	    }
	} else {
	    Rsteqr("I", matsiz, &d[submat], &e[submat], &work[iq - 1 + iwork[iqptr + curr]], matsiz, &work[0], info);
	    if (*info != 0) {
		goto L130;
	    }
	    if (icompq == 1) {
		Rgemm("N", "N", qsiz, matsiz, matsiz, One, &q[submat * ldq + 1], ldq, &work[iq - 1 + iwork[iqptr + curr]], matsiz, Zero, &qstore[submat * ldqs + 1], ldqs);
	    }
	    iwork[iqptr + curr + 1] = iwork[iqptr + curr] + matsiz * matsiz;
	    ++curr;
	}
	k = 0;
	for (j = submat; j <= iwork[i + 1]; j++) {
	    iwork[indxq + j] = k;
	    k++;
	}
    }
//Successively merge eigensystems of adjacent submatrices
//into eigensystem for the corresponding larger matrix.
//while ( SUBPBS > 1 )
    curlvl = 0;

  L80:
    if (subpbs > 1) {
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
//into an eigensystem of size MATSIZ.
//DLAED1 is used only for the full eigensystem of a tridiagonal
//matrix.
//DLAED7 handles the cases in which eigenvalues only or eigenvalues
//and eigenvectors of a full symmetric matrix (which was reduced to
//tridiagonal form) are desired.
	    if (icompq == 2) {
		Rlaed1(matsiz, &d[submat], &q[submat + submat * ldq], ldq, &iwork[indxq + submat], e[submat + msd2 - 1], msd2, &work[0], &iwork[subpbs + 1], info);
	    } else {
		Rlaed7(icompq, matsiz, qsiz, tlvls, curlvl, curprb, &d[submat],
		       &qstore[submat * ldqs + 1], ldqs, &iwork[indxq + submat], e[submat + msd2 - 1], msd2, &work[iq], &iwork[iqptr], &iwork[iprmpt], &iwork[iperm]
		       , &iwork[igivpt], &iwork[igivcl], &work[igivnm], &work[iwrem], &iwork[subpbs + 1], info);
	    }
	    if (*info != 0) {
		goto L130;
	    }
	    iwork[i / 2 + 1] = iwork[i + 2];

	}
	subpbs /= 2;
	++curlvl;
	goto L80;
    }
//end while
//Re-merge the eigenvalues/vectors which were deflated at the final
//merge step.
    if (icompq == 1) {
	for (i = 0; i < n; i++) {
	    j = iwork[indxq + i];
	    work[i] = d[j];
	    Rcopy(qsiz, &qstore[j * ldqs + 1], 1, &q[i * ldq + 1], 1);

	}
	Rcopy(n, &work[0], 1, &d[0], 1);
    } else if (icompq == 2) {
	for (i = 0; i < n; i++) {
	    j = iwork[indxq + i];
	    work[i] = d[j];
	    Rcopy(n, &q[j * ldq + 1], 1, &work[n * i + 1], 1);
	}
	Rcopy(n, &work[0], 1, &d[0], 1);
	Rlacpy("A", n, n, &work[n + 1], n, &q[0], ldq);
    } else {
	for (i = 0; i < n; i++) {
	    j = iwork[indxq + i];
	    work[i] = d[j];
	}
	Rcopy(n, &work[0], 1, &d[0], 1);
    }
    goto L140;

  L130:
    *info = submat * (n + 1) + submat + matsiz - 1;
  L140:
    return;
}
