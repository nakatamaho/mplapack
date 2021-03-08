/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd3.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasd3(INTEGER nl, INTEGER nr, INTEGER sqre, INTEGER k, REAL * d, REAL * q, INTEGER ldq,
       REAL * dsigma, REAL * u, INTEGER ldu, REAL * u2, INTEGER ldu2, REAL * vt, INTEGER ldvt, REAL * vt2, INTEGER ldvt2, INTEGER * idxc, INTEGER * ctot, REAL * z, INTEGER * info)
{
    INTEGER i, j, m, n, jc;
    REAL rho;
    INTEGER nlp1, nlp2, nrp1;
    REAL temp;
    INTEGER ctemp;
    INTEGER ktemp;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
    *info = 0;
    if (nl < 1) {
	*info = -1;
    } else if (nr < 1) {
	*info = -2;
    } else if (sqre != 1 && sqre != 0) {
	*info = -3;
    }

    n = nl + nr + 1;
    m = n + sqre;
    nlp1 = nl + 1;
    nlp2 = nl + 2;

    if (k < 1 || k > n) {
	*info = -4;
    } else if (ldq < k) {
	*info = -7;
    } else if (ldu < n) {
	*info = -10;
    } else if (ldu2 < n) {
	*info = -12;
    } else if (ldvt < m) {
	*info = -14;
    } else if (ldvt2 < m) {
	*info = -16;
    }
    if (*info != 0) {
	Mxerbla("Rlasd3", -(*info));
	return;
    }
//Quick return if possible
    if (k == 1) {
	d[1] = abs(z[1]);
	Rcopy(m, &vt2[ldvt2 + 1], ldvt2, &vt[ldvt + 1], ldvt);
	if (z[1] > Zero) {
	    Rcopy(n, &u2[ldu2 + 1], 1, &u[ldu + 1], 1);
	} else {
	    for (i = 0; i < n; i++) {
		u[i + ldu] = -u2[i + ldu2];
	    }
	}
	return;
    }
//Keep a copy of Z.
    Rcopy(k, &z[1], 1, &q[0], 1);

//Normalize Z.
    rho = Rnrm2(k, &z[1], 1);
    Rlascl("G", 0, 0, rho, One, k, 1, &z[1], k, info);
    rho *= rho;
//Find the new singular values.
    for (j = 0; j < k; j++) {
	Rlasd4(k, j, &dsigma[1], &z[1], &u[j * ldu + 1], rho, &d[j], &vt[j * ldvt + 1], info);
//If the zero finder fails, the computation is terminated.
	if (*info != 0) {
	    return;
	}
    }
//Compute updated Z.
    for (i = 0; i < k; i++) {
	z[i] = u[i + k * ldu] * vt[i + k * ldvt];
	for (j = 0; j < i - 1; j++) {
	    z[i] *= u[i + j * ldu] * vt[i + j * ldvt] / (dsigma[i] - dsigma[j]) / (dsigma[i] + dsigma[j]);
	}
	for (j = i; j <= k - 1; j++) {
	    z[i] *= u[i + j * ldu] * vt[i + j * ldvt] / (dsigma[i] - dsigma[j + 1]) / (dsigma[i] + dsigma[j + 1]);
	}
	z[i] = sign(sqrt(abs(z[i])), q[i + ldq]);
    }
//Compute left singular vectors of the modified diagonal matrix,
//and store related information for the right singular vectors.
    for (i = 0; i < k; i++) {
	vt[i * ldvt + 1] = z[1] / u[i * ldu + 1] / vt[i * ldvt + 1];
	u[i * ldu + 1] = -One;
	for (j = 2; j <= k; j++) {
	    vt[j + i * ldvt] = z[j] / u[j + i * ldu] / vt[j + i * ldvt];
	    u[j + i * ldu] = dsigma[j] * vt[j + i * ldvt];
	}
	temp = Rnrm2(k, &u[i * ldu + 1], 1);
	q[i * ldq + 1] = u[i * ldu + 1] / temp;
	for (j = 2; j <= k; j++) {
	    jc = idxc[j];
	    q[j + i * ldq] = u[jc + i * ldu] / temp;
	}
    }
//Update the left singular vector matrix.
    if (k == 2) {
	Rgemm("N", "N", n, k, k, One, &u2[0], ldu2, &q[0], ldq, Zero, &u[0], ldu);
	goto L100;
    }
    if (ctot[1] > 0) {
	Rgemm("N", "N", nl, k, ctot[1], One, &u2[(ldu2 << 1) + 1], ldu2, &q[ldq + 2], ldq, Zero, &u[ldu + 1], ldu);
	if (ctot[3] > 0) {
	    ktemp = ctot[1] + 2 + ctot[2];
	    Rgemm("N", "N", nl, k, ctot[3], One, &u2[ktemp * ldu2 + 1]
		  , ldu2, &q[ktemp + ldq], ldq, One, &u[ldu + 1], ldu);
	}
    } else if (ctot[3] > 0) {
	ktemp = ctot[1] + 2 + ctot[2];
	Rgemm("N", "N", nl, k, ctot[3], One, &u2[ktemp * ldu2 + 1], ldu2, &q[ktemp + ldq], ldq, Zero, &u[ldu + 1], ldu);
    } else {
	Rlacpy("F", nl, k, &u2[0], ldu2, &u[0], ldu);
    }
    Rcopy(k, &q[ldq + 1], ldq, &u[nlp1 + ldu], ldu);
    ktemp = ctot[1] + 2;
    ctemp = ctot[2] + ctot[3];
    Rgemm("N", "N", nr, k, ctemp, One, &u2[nlp2 + ktemp * ldu2], ldu2, &q[ktemp + ldq], ldq, Zero, &u[nlp2 + ldu], ldu);
//Generate the right singular vectors.
  L100:
    for (i = 0; i < k; i++) {
	temp = Rnrm2(k, &vt[i * ldvt + 1], 1);
	q[i + ldq] = vt[i * ldvt + 1] / temp;
	for (j = 2; j <= k; j++) {
	    jc = idxc[j];
	    q[i + j * ldq] = vt[jc + i * ldvt] / temp;
	}
    }
//Update the right singular vector matrix.
    if (k == 2) {
	Rgemm("N", "N", k, m, k, One, &q[0], ldq, &vt2[0], ldvt2, Zero, &vt[0], ldvt);
	return;
    }
    ktemp = ctot[1] + 1;
    Rgemm("N", "N", k, nlp1, ktemp, One, &q[ldq + 1], ldq, &vt2[ldvt2 + 1], ldvt2, Zero, &vt[ldvt + 1], ldvt);
    ktemp = ctot[1] + 2 + ctot[2];
    if (ktemp <= ldvt2) {
	Rgemm("N", "N", k, nlp1, ctot[3], One, &q[ktemp * ldq + 1], ldq, &vt2[ktemp + ldvt2], ldvt2, One, &vt[ldvt + 1], ldvt);
    }

    ktemp = ctot[1] + 1;
    nrp1 = nr + sqre;
    if (ktemp > 1) {
	for (i = 0; i < k; i++) {
	    q[i + ktemp * ldq] = q[i + ldq];

	}
	for (i = nlp2; i <= m; i++) {
	    vt2[ktemp + i * ldvt2] = vt2[i * ldvt2 + 1];
	}
    }
    ctemp = ctot[2] + 1 + ctot[3];
    Rgemm("N", "N", k, nrp1, ctemp, One, &q[ktemp * ldq + 1], ldq, &vt2[ktemp + nlp2 * ldvt2], ldvt2, Zero, &vt[nlp2 * ldvt + 1], ldvt);
    return;
}
