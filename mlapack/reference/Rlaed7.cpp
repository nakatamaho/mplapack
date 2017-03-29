/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed7.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaed7(INTEGER icompq, INTEGER n, INTEGER qsiz, INTEGER tlvls, INTEGER curlvl, INTEGER curpbm,
       REAL * d, REAL * q, INTEGER ldq, INTEGER * indxq, REAL rho,
       INTEGER cutpnt, REAL * qstore, INTEGER * qptr, INTEGER * prmptr, INTEGER * perm,
       INTEGER * givptr, INTEGER * givcol, REAL * givnum, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, k, is, iw, iz, iq2, ptr, ldq2, indx, curr;
    INTEGER indxc, indxp;
    INTEGER idlmda;
    INTEGER coltyp;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (icompq == 1 && qsiz < n) {
	*info = -4;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -9;
    } else if (max((INTEGER) 1, n) > cutpnt || n < cutpnt) {
	*info = -12;
    }
    if (*info != 0) {
	Mxerbla("Rlaed7", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//The following values are for bookkeeping purposes only.  They are
//INTEGER pointers which indicate the portion of the workspace
//used by a particular array in DLAED8 and DLAED9.
    if (icompq == 1) {
	ldq2 = qsiz;
    } else {
	ldq2 = n;
    }

    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq2 = iw + n;
    is = iq2 + n * ldq2;
    indx = 1;
    indxc = indx + n;
    coltyp = indxc + n;
    indxp = coltyp + n;
//Form the z-vector which consists of the last row of Q_1 and the
//first row of Q_2
    ptr = 1 + (2 ^ tlvls);
    for (i = 0; i < curlvl - 1; i++) {
	ptr = ptr + (2 ^ (tlvls - i));
    }
    curr = ptr + curpbm;
    Rlaeda(n, tlvls, curlvl, curpbm, &prmptr[1], &perm[1], &givptr[1], &givcol[3], &givnum[3], &qstore[1], &qptr[1], &work[iz], &work[iz + n], info);

//When solving the final problem, we no longer need the stored data,
//so we will overwrite the data from this level onto the previously
//used storage space.
    if (curlvl == tlvls) {
	qptr[curr] = 1;
	prmptr[curr] = 1;
	givptr[curr] = 1;
    }
//Sort and Deflate eigenvalues.
    Rlaed8(icompq, &k, n, qsiz, &d[0], &q[0], ldq, &indxq[1], &rho,
	   cutpnt, &work[iz], &work[idlmda], &work[iq2], ldq2, &work[iw], &perm[prmptr[curr]], &givptr[curr + 1], &givcol[(givptr[curr] << 1)
															  + 1], &givnum[(givptr[curr] << 1) + 1], &iwork[indxp],
	   &iwork[indx], info);
    prmptr[curr + 1] = prmptr[curr] + n;
    givptr[curr + 1] += givptr[curr];
//Solve Secular Equation.
    if (k != 0) {
	Rlaed9(k, 1, k, n, &d[0], &work[is], k, rho, &work[idlmda], &work[iw], &qstore[qptr[curr]], k, info);
	if (*info != 0) {
	    goto L30;
	}
	if (icompq == 1) {
	    Rgemm("N", "N", qsiz, k, k, One, &work[iq2], ldq2, &qstore[qptr[curr]], k, Zero, &q[0], ldq);
	}
	qptr[curr + 1] = qptr[curr] + k * k;
//Prepare the INDXQ sorting permutation.
	Rlamrg(k, n - k, &d[0], 1, -1, &indxq[1]);
    } else {
	qptr[curr + 1] = qptr[curr];
	for (i = 0; i < n; i++) {
	    indxq[i] = i;
	}
    }
  L30:
    return;
}
