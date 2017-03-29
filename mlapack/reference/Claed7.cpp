/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claed7.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Claed7(INTEGER n, INTEGER cutpnt, INTEGER qsiz, INTEGER tlvls, INTEGER curlvl, INTEGER curpbm,
       REAL * d, COMPLEX * q, INTEGER ldq, REAL rho, INTEGER * indxq,
       REAL * qstore, INTEGER * qptr, INTEGER * prmptr, INTEGER * perm, INTEGER * givptr,
       INTEGER * givcol, REAL * givnum, COMPLEX * work, REAL * rwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, k, iq, iw, iz, ptr, indx, curr, indxc, indxp, idlmda;
    INTEGER coltyp;

//Test the input parameters.
    *info = 0;
    if (n < 0) {
	*info = -1;
    } else if (max((INTEGER) 1, n) > cutpnt || n < cutpnt) {
	*info = -2;
    } else if (qsiz < n) {
	*info = -3;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -9;
    }
    if (*info != 0) {
	Mxerbla("Claed7", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//The following values are for bookkeeping purposes only.  They are
//INTEGER pointers which indicate the portion of the workspace
//used by a particular array in DLAED2 and SLAED3.
    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq = iw + n;
    indx = 1;
    indxc = indx + n;
    coltyp = indxc + n;
    indxp = coltyp + n;
//Form the z-vector which consists of the last row of Q_1 and the
//first row of Q_2
    ptr = (2 ^ tlvls) + 1;
    for (i = 0; i < curlvl - 1; i++) {
	ptr += (2 ^ (tlvls - i));
    }
    curr = ptr + curpbm;
    Rlaeda(n, tlvls, curlvl, curpbm, &prmptr[1], &perm[1], &givptr[1], &givcol[3], &givnum[3], &qstore[1], &qptr[1], &rwork[iz], &rwork[iz + n], info);
//When solving the final problem, we no longer need the stored data,
//so we will overwrite the data from this level onto the previously
//used storage space.
    if (curlvl == tlvls) {
	qptr[curr] = 1;
	prmptr[curr] = 1;
	givptr[curr] = 1;
    }
//Sort and Deflate eigenvalues.
    Claed8(&k, n, qsiz, &q[0], ldq, &d[0], &rho, cutpnt, &rwork[iz],
	   &rwork[idlmda], &work[0], qsiz, &rwork[iw], &iwork[indxp],
	   &iwork[indx], &indxq[1], &perm[prmptr[curr]], &givptr[curr + 1], &givcol[(givptr[curr] << 1) + 1], &givnum[(givptr[curr] << 1) + 1], info);
    prmptr[curr + 1] = prmptr[curr] + n;
    givptr[curr + 1] += givptr[curr];
//Solve Secular Equation.
    if (k != 0) {
	Rlaed9(k, 1, k, n, &d[0], &rwork[iq], k, rho, &rwork[idlmda]
	       , &rwork[iw], &qstore[qptr[curr]], k, info);
	Clacrm(qsiz, k, &work[0], qsiz, &qstore[qptr[curr]], k, &q[0], ldq, &rwork[iq]);
	qptr[curr + 1] = qptr[curr] + k * k;
	if (*info != 0) {
	    return;
	}
//Prepare the INDXQ sorting premutation.
	Rlamrg(k, n - k, &d[0], 1, -1, &indxq[1]);
    } else {
	qptr[curr + 1] = qptr[curr];
	for (i = 0; i < n; i++) {
	    indxq[i] = i;
	}
    }
    return;
}
