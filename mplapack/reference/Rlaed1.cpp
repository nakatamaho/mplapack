/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed1.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlaed1(INTEGER n, REAL * d, REAL * q, INTEGER ldq, INTEGER * indxq, REAL rho, INTEGER cutpnt, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, k, is, iw, iz, iq2, zpp1, indx, indxc;
    INTEGER indxp;
    INTEGER idlmda;
    INTEGER coltyp;

//Test the input parameters.
    *info = 0;
    if (n < 0) {
	*info = -1;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -4;
    } else {
	if (max((INTEGER) 1, n / 2) > cutpnt || n / 2 < cutpnt) {
	    *info = -7;
	}
    }
    if (*info != 0) {
	Mxerbla("Rlaed1", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

//The following values are INTEGER pointers which indicate
//the portion of the workspace
//used by a particular array in DLAED2 and DLAED3.
    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq2 = iw + n;

    indx = 1;
    indxc = indx + n;
    coltyp = indxc + n;
    indxp = coltyp + n;

//Form the z-vector which consists of the last row of Q_1 and the
//first row of Q_Two

    Rcopy(cutpnt, &q[cutpnt + ldq], ldq, &work[iz], 1);
    zpp1 = cutpnt + 1;
    Rcopy(n - cutpnt, &q[zpp1 + zpp1 * ldq], ldq, &work[iz + cutpnt], 1);

//Deflate eigenvalues.
    Rlaed2(&k, n, cutpnt, &d[0], &q[0], ldq, &indxq[1], &rho, &work[iz], &work[idlmda], &work[iw], &work[iq2], &iwork[indx], &iwork[indxc], &iwork[indxp], &iwork[coltyp], info);

    if (*info != 0)
	return;

//Solve Secular Equation.

    if (k != 0) {
	is = (iwork[coltyp] + iwork[coltyp + 1]) * cutpnt + (iwork[coltyp + 1] + iwork[coltyp + 2]) * (n - cutpnt) + iq2;
	Rlaed3(k, n, cutpnt, &d[0], &q[0], ldq, &rho, &work[idlmda], &work[iq2], &iwork[indxc], &iwork[coltyp], &work[iw], &work[is], info);
	if (*info != 0)
	    return;

//Prepare the INDXQ sorting permutation.
	Rlamrg(k, n - k, &d[0], 1, -1, &indxq[1]);
    } else {
	for (i = 0; i < n; i++) {
	    indxq[i] = i;

	}
    }
    return;
}
