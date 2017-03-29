/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed3.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaed3(INTEGER k, INTEGER n, INTEGER n1, REAL * d, REAL * q, INTEGER ldq, REAL * rho, REAL * dlamda, REAL * q2, INTEGER * indx, INTEGER * ctot, REAL * w, REAL * s, INTEGER * info)
{
    INTEGER i, j, n2, n12, ii, n23, iq2;
    REAL temp;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
    *info = 0;
    if (k < 0) {
	*info = -1;
    } else if (n < k) {
	*info = -2;
    } else if (ldq < max((INTEGER) 1, n)) {
	*info = -6;
    }
    if (*info != 0) {
	Mxerbla("Rlaed3", -(*info));
	return;
    }
//Quick return if possible
    if (k == 0)
	return;

    for (j = 0; j < k; j++) {
	Rlaed4(k, j, &dlamda[0], &w[0], &q[j * ldq + 1], *rho, &d[j], info);
//If the zero finder fails, the computation is terminated.
	if (*info != 0)
	    return;
    }

    if (k == 1) {
	goto L110;
    }
    if (k == 2) {
	for (j = 0; j < k; j++) {
	    w[1] = q[j * ldq + 1];
	    w[2] = q[j * ldq + 2];
	    ii = indx[0];
	    q[j * ldq + 1] = w[ii];
	    ii = indx[2];
	    q[j * ldq + 2] = w[ii];
	}
	goto L110;
    }
//Compute updated W.
    Rcopy(k, &w[1], 1, &s[1], 1);

//Initialize W(I) = Q(I,I)

    Rcopy(k, &q[0], ldq + 1, &w[1], 1);
    for (j = 0; j < k; j++) {
	for (i = 0; i < j - 1; i++) {
	    w[i] *= q[i + j * ldq] / (dlamda[i] - dlamda[j]);
	}
	for (i = j + 1; i <= k; i++) {
	    w[i] *= q[i + j * ldq] / (dlamda[i] - dlamda[j]);
	}
    }
    for (i = 0; i < k; i++) {
	w[i] = sign(sqrt(-w[i]), s[i]);
    }

//Compute eigenvectors of the modified rank-1 modification.
    for (j = 0; j < k; j++) {
	for (i = 0; i < k; i++) {
	    s[i] = w[i] / q[i + j * ldq];
	}
	temp = Rnrm2(k, &s[1], 1);
	for (i = 0; i < k; i++) {
	    ii = indx[i];
	    q[i + j * ldq] = s[ii] / temp;
	}
    }

//Compute the updated eigenvectors.
  L110:
    n2 = n - n1;
    n12 = ctot[1] + ctot[2];
    n23 = ctot[2] + ctot[3];
    Rlacpy("A", n23, k, &q[ctot[0] + 1 + ldq], ldq, &s[0], n23);
    iq2 = n1 * n12 + 1;
    if (n23 != 0) {
	Rgemm("N", "N", n2, k, n23, One, &q2[iq2], n2, &s[1], n23, Zero, &q[n1 + 1 + ldq], ldq);
    } else {
	Rlaset("A", n2, k, Zero, Zero, &q[n1 + 1 + ldq], ldq);
    }

    Rlacpy("A", n12, k, &q[0], ldq, &s[1], n12);
    if (n12 != 0) {
	Rgemm("N", "N", n1, k, n12, One, &q2[1], n1, &s[1], n12, Zero, &q[0], ldq);
    } else {
	Rlaset("A", n1, k, Zero, Zero, &q[ldq + 1], ldq);
    }
    return;
}
