/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd1.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasd1(INTEGER nl, INTEGER nr, INTEGER * sqre, REAL * d, REAL * alpha,
       REAL * beta, REAL * u, INTEGER ldu, REAL * vt, INTEGER ldvt, INTEGER * idxq, INTEGER * iwork, REAL * work, INTEGER * info)
{
    INTEGER i, k, m, n, iq, iz, iu2, ldq, idx, ldu2, ivt2, idxc, idxp, ldvt2;
    INTEGER isigma;
    REAL orgnrm;
    INTEGER coltyp;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    if (nl < 1) {
	*info = -1;
    } else if (nr < 1) {
	*info = -2;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -3;
    }
    if (*info != 0) {
	Mxerbla("Rlasd1", -(*info));
	return;
    }
    n = nl + nr + 1;
    m = n + *sqre;

//The following values are for bookkeeping purposes only.  They are
//INTEGER pointers which indicate the portion of the workspace
//used by a particular array in DLASD2 and DLASD3.
    ldu2 = n;
    ldvt2 = m;

    iz = 1;
    isigma = iz + m;
    iu2 = isigma + n;
    ivt2 = iu2 + ldu2 * n;
    iq = ivt2 + ldvt2 * m;

    idx = 1;
    idxc = idx + n;
    coltyp = idxc + n;
    idxp = coltyp + n;

//Scale.
    mtemp1 = abs(*alpha), mtemp2 = abs(*beta);
    orgnrm = max(mtemp1, mtemp2);
    d[nl + 1] = Zero;
    for (i = 0; i < n; i++) {
	if (abs(d[i]) > orgnrm) {
	    orgnrm = abs(d[i]);
	}
    }
    Rlascl("G", 0, 0, orgnrm, One, n, 1, &d[0], n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

//Deflate singular values.
    Rlasd2(nl, nr, *sqre, &k, &d[0], &work[iz], *alpha, *beta, &u[0],
	   ldu, &vt[0], ldvt, &work[isigma], &work[iu2], ldu2, &work[ivt2], ldvt2, &iwork[idxp], &iwork[idx], &iwork[idxc], &idxq[1], &iwork[coltyp], info);
//Solve Secular Equation and update singular vectors.
    ldq = k;
    Rlasd3(nl, nr, *sqre, k, &d[0], &work[iq], ldq, &work[isigma], &u[0], ldu, &work[iu2], ldu2, &vt[0], ldvt, &work[ivt2], ldvt2, &iwork[idxc], &iwork[coltyp], &work[iz], info);
    if (*info != 0) {
	return;
    }
//Unscale.
    Rlascl("G", 0, 0, One, orgnrm, n, 1, &d[0], n, info);

//Prepare the IDXQ sorting permutation.
    Rlamrg(k, n - k, &d[0], 1, -1, &idxq[1]);
    return;
}
