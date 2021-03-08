/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd6.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasd6(INTEGER icompq, INTEGER nl, INTEGER nr, INTEGER sqre, REAL * d, REAL * vf,
       REAL * vl, REAL * alpha, REAL * beta, INTEGER * idxq, INTEGER * perm,
       INTEGER * givptr, INTEGER * givcol, INTEGER ldgcol, REAL * givnum, INTEGER ldgnum,
       REAL * poles, REAL * difl, REAL * difr, REAL * z, INTEGER k, REAL * c, REAL * s, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, m, n, iw, idx, idxc, idxp, ivfw, ivlw;
    INTEGER isigma;
    REAL orgnrm;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    n = nl + nr + 1;
    m = n + sqre;

    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (nl < 1) {
	*info = -2;
    } else if (nr < 1) {
	*info = -3;
    } else if (sqre < 0 || sqre > 1) {
	*info = -4;
    } else if (ldgcol < n) {
	*info = -14;
    } else if (ldgnum < n) {
	*info = -16;
    }
    if (*info != 0) {
	Mxerbla("Rlasd6", -(*info));
	return;
    }
//The following values are for bookkeeping purposes only.  They are
//INTEGER pointers which indicate the portion of the workspace
//used by a particular array in DLASD7 and DLASD8.
    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;

    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;
//Scale.
    orgnrm = max(abs(*alpha), abs(*beta));
    d[nl + 1] = Zero;
    for (i = 0; i < n; i++) {
	if (abs(d[i]) > orgnrm) {
	    orgnrm = abs(d[i]);
	}

    }
    Rlascl("G", 0, 0, orgnrm, One, n, 1, &d[0], n, info);
    *alpha = *alpha / orgnrm;
    *beta = *beta / orgnrm;
//Sort and Deflate singular values.
    Rlasd7(icompq, nl, nr, sqre, k, &d[0], &z[1], &work[iw], &vf[1],
	   &work[ivfw], &vl[1], &work[ivlw], *alpha, *beta, &work[isigma],
	   &iwork[idx], &iwork[idxp], &idxq[1], &perm[1], givptr, &givcol[0], ldgcol, &givnum[0], ldgnum, c, s, info);
//Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.
    Rlasd8(icompq, k, &d[0], &z[1], &vf[1], &vl[1], &difl[1], &difr[1], ldgnum, &work[isigma], &work[iw], info);
//Save the poles if ICOMPQ = One
    if (icompq == 1) {
	Rcopy(k, &d[0], 1, &poles[ldgnum + 1], 1);
	Rcopy(k, &work[isigma], 1, &poles[(ldgnum * 2) + 1], 1);
    }
//Unscale.
    Rlascl("G", 0, 0, One, orgnrm, n, 1, &d[0], n, info);
//Prepare the IDXQ sorting permutation.
    Rlamrg(k, n - k, &d[0], 1, -1, &idxq[1]);
    return;
}
