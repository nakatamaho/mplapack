/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cstein.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Cstein(INTEGER n, REAL * d, REAL * e, INTEGER m, REAL * w, INTEGER * iblock, INTEGER * isplit, COMPLEX * z, INTEGER ldz,
       REAL * work, INTEGER * iwork, INTEGER * ifail, INTEGER * info)
{
    INTEGER i, j, b1, j1, bn, jr;
    REAL xj, scl, eps, sep, nrm, tol;
    INTEGER its;
    REAL xjm = 0, ztr, eps1;
    INTEGER jblk, nblk, jmax;
    INTEGER iseed[4], gpind = 0, iinfo;
    REAL ortol = 0;
    INTEGER indrv1, indrv2, indrv3, indrv4, indrv5;
    INTEGER nrmchk;
    INTEGER blksiz;
    REAL onenrm = 0.0, dtpcrt = 0.0, pertol;
    REAL Zero = 0.0, One = 1.0, Ten = 10.0, ODM3 = 0.001, ODM1 = 0.1;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    for (i = 0; i < m; i++) {
	ifail[i] = 0;
    }
    if (n < 0) {
	*info = -1;
    } else if (m < 0 || m > n) {
	*info = -4;
    } else if (ldz < max((INTEGER) 1, n)) {
	*info = -9;
    } else {
	for (j = 2; j <= m; j++) {
	    if (iblock[j] < iblock[j - 1]) {
		*info = -6;
		goto L30;
	    }
	    if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
		*info = -5;
		goto L30;
	    }

	}
      L30:
	;
    }
    if (*info != 0) {
	Mxerbla("Cstein", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0 || m == 0) {
	return;
    } else if (n == 1) {
	z[ldz + 1] = One;
	return;
    }
//Get machine constants.
    eps = Rlamch("Precision");
//Initialize seed for random number generator DLARNV.
    for (i = 1; i < 4; i++) {
	iseed[i - 1] = 1;
    }
//Initialize poINTEGERers.
    indrv1 = 0;
    indrv2 = indrv1 + n;
    indrv3 = indrv2 + n;
    indrv4 = indrv3 + n;
    indrv5 = indrv4 + n;
//Compute eigenvectors of matrix blocks.
    j1 = 1;
    for (nblk = 0; nblk <= iblock[m]; ++nblk) {
//Find starting and ending indices of block nblk.
	if (nblk == 1) {
	    b1 = 1;
	} else {
	    b1 = isplit[nblk - 1] + 1;
	}
	bn = isplit[nblk];
	blksiz = bn - b1 + 1;
	if (blksiz == 1) {
	    goto L60;
	}
	gpind = b1;
//Compute reorthogonalization criterion and stopping criterion.
	onenrm = abs(d[b1]) + abs(e[b1]);
	mtemp1 = onenrm, mtemp2 = abs(d[bn]) + abs(e[bn - 1]);
	onenrm = max(mtemp1, mtemp2);
	for (i = b1 + 1; i <= bn - 1; i++) {
	    mtemp1 = onenrm, mtemp2 = abs(d[i]) + abs(e[i - 1]) + abs(e[i]);
	    onenrm = max(mtemp1, mtemp2);
	}
	ortol = onenrm * ODM3;
	dtpcrt = sqrt(ODM1 / blksiz);
//Loop through eigenvalues of block nblk.
      L60:
	jblk = 0;
	for (j = j1; j <= m; j++) {
	    if (iblock[j] != nblk) {
		j1 = j;
		goto L180;
	    }
	    jblk++;
	    xj = w[j];
//Skip all the work if the block size is one.
	    if (blksiz == 1) {
		work[indrv1 + 1] = One;
		goto L140;
	    }
//If eigenvalues j and j-1 are too close, add a relatively
//small perturbation.
	    if (jblk > 1) {
		eps1 = abs(eps * xj);
		pertol = eps1 * Ten;
		sep = xj - xjm;
		if (sep < pertol) {
		    xj = xjm + pertol;
		}
	    }
	    its = 0;
	    nrmchk = 0;
//Get random starting vector.
	    Rlarnv(2, iseed, blksiz, &work[indrv1 + 1]);
//Copy the matrix T so it won't be destroyed in factorization.
	    Rcopy(blksiz, &d[b1], 1, &work[indrv4 + 1], 1);
	    Rcopy(blksiz - 1, &e[b1], 1, &work[indrv2 + 2], 1);
	    Rcopy(blksiz - 1, &e[b1], 1, &work[indrv3 + 1], 1);
//Compute LU factors with partial pivoting  ( PT = LU )
	    tol = Zero;
	    Rlagtf(blksiz, &work[indrv4 + 1], xj, &work[indrv2 + 2], &work[indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);
//Update iteration count.
	  L70:
	    its++;
	    if (its > 5) {
		goto L120;
	    }
//Normalize and scale the righthand side vector Pb.
	    mtemp1 = eps, mtemp2 = abs(work[indrv4 + blksiz]);
	    scl = blksiz * onenrm * max(mtemp1, mtemp2) / Rasum(blksiz, &work[indrv1 + 1], 1);
	    Rscal(blksiz, scl, &work[indrv1 + 1], 1);
//Solve the system LU = Pb.
	    Rlagts(-1, blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[indrv1 + 1], &tol, &iinfo);
//Reorthogonalize by modified Gram-Schmidt if eigenvalues are
//close enough.
	    if (jblk == 1) {
		goto L110;
	    }
	    if (abs(xj - xjm) > ortol) {
		gpind = j;
	    }
	    if (gpind != j) {
		for (i = gpind; i <= j - 1; i++) {
		    ztr = Zero;
		    for (jr = 1; jr <= blksiz; jr++) {
			ztr = ztr + work[indrv1 + jr] * z[b1 - 1 + jr + i * ldz].real();

		    }
		    for (jr = 1; jr <= blksiz; jr++) {
			work[indrv1 + jr] = work[indrv1 + jr] - ztr * z[b1 - 1 + jr + i * ldz].real();
		    }
		}
	    }
//Check the infinity norm of the iterate.
	  L110:
	    jmax = iRamax(blksiz, &work[indrv1 + 1], 1);
	    nrm = abs(work[indrv1 + jmax]);
//Continue for additional iterations after norm reaches
//stopping criterion.
	    if (nrm < dtpcrt) {
		goto L70;
	    }
	    ++nrmchk;
	    if (nrmchk < 3) {
		goto L70;
	    }
	    goto L130;
//If stopping criterion was not satisfied, update info and
//store eigenvector number in array ifail.
	  L120:
	    ++(*info);
	    ifail[*info] = j;
//Accept iterate as jth eigenvector.
	  L130:
	    scl = Zero / Rnrm2(blksiz, &work[indrv1 + 1], 1);
	    jmax = iRamax(blksiz, &work[indrv1 + 1], 1);
	    if (work[indrv1 + jmax] < Zero) {
		scl = -scl;
	    }
	    Rscal(blksiz, scl, &work[indrv1 + 1], 1);
	  L140:
	    for (i = 0; i < n; i++) {
		z[i + j * ldz] = Zero;

	    }
	    for (i = 0; i < blksiz; i++) {
		z[b1 + i - 1 + j * ldz] = (COMPLEX) work[indrv1 + i];
	    }
//Save the shift to check eigenvalue spacing at next
//iteration.
	    xjm = xj;
	}
      L180:
	;
    }
    return;
}
