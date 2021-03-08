/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctrsen.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctrsen(const char *job, const char *compq, LOGICAL * select, INTEGER
	    n, COMPLEX * t, INTEGER ldt, COMPLEX * q, INTEGER ldq, COMPLEX * w, INTEGER m, REAL * s, REAL * sep, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER k, n1, n2, nn, ks;
    REAL est;
    INTEGER kase, ierr;
    REAL scale;
    INTEGER isave[3], lwmin = 0;
    INTEGER wantq, wants;
    REAL rnorm, rwork[1];
    INTEGER wantbh;
    INTEGER wantsp;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;

//Decode and test the input parameters.
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantsp = Mlsame(job, "V") || wantbh;
    wantq = Mlsame(compq, "V");
//Set M to the number of selected eigenvalues.
    m = 0;
    for (k = 0; k < n; k++) {
	if (select[k]) {
	    ++(m);
	}
    }

    n1 = m;
    n2 = n - m;
    nn = n1 * n2;
    *info = 0;
    lquery = lwork == -1;

    if (wantsp) {
	lwmin = max((INTEGER) 1, nn * 2);
    } else if (Mlsame(job, "N")) {
	lwmin = 1;
    } else if (Mlsame(job, "E")) {
	lwmin = max((INTEGER) 1, nn);
    }
    if (!Mlsame(job, "N") && !wants && !wantsp) {
	*info = -1;
    } else if (!Mlsame(compq, "N") && !wantq) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (ldt < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldq < 1 || (wantq && ldq < n)) {
	*info = -8;
    } else if (lwork < lwmin && !lquery) {
	*info = -14;
    }
    if (*info == 0) {
	work[1] = lwmin;
    }
    if (*info != 0) {
	Mxerbla("Ctrsen", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == n || m == 0) {
	if (wants) {
	    *s = One;
	}
	if (wantsp) {
	    *sep = Clange("1", n, n, &t[0], ldt, rwork);
	}
	goto L40;
    }
//Collect the selected eigenvalues at the top left corner of T.
    ks = 0;
    for (k = 0; k < n; k++) {
	if (select[k]) {
	    ks++;
//Swap the K-th eigenvalue to position KS.
	    if (k != ks) {
		Ctrexc(compq, n, &t[0], ldt, &q[0], ldq, k, ks, &ierr);
	    }
	}
    }
    if (wants) {
//Solve the Sylvester equation for R:
//   T11*R - R*T22 = scale*T12
	Clacpy("F", n1, n2, &t[(n1 + 1) * ldt + 1], ldt, &work[0], n1);
	Ctrsyl("N", "N", -1, n1, n2, &t[0], ldt, &t[n1 + 1 + (n1 + 1) * ldt], ldt, &work[0], n1, &scale, &ierr);
//Estimate the reciprocal of the condition number of the cluster
//of eigenvalues.
	rnorm = Clange("F", n1, n2, &work[0], n1, rwork);
	if (rnorm == Zero) {
	    *s = One;
	} else {
	    *s = scale / (sqrt(scale * scale / rnorm + rnorm) * sqrt(rnorm));
	}
    }
    if (wantsp) {
//Estimate sep(T11,T22).
	est = Zero;
	kase = 0;
      L30:
	Clacn2(nn, &work[nn + 1], &work[0], &est, &kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
//Solve T11*R - R*T22 = scale*X.
		Ctrsyl("N", "N", -1, n1, n2, &t[0], ldt, &t[n1 + 1 + (n1 + 1) * ldt], ldt, &work[0], n1, &scale, &ierr);
	    } else {
//Solve T11'*R - R*T22' = scale*X.
		Ctrsyl("C", "C", -1, n1, n2, &t[0], ldt, &t[n1 + 1 + (n1 + 1) * ldt], ldt, &work[0], n1, &scale, &ierr);
	    }
	    goto L30;
	}
	*sep = scale / est;
    }
  L40:
//Copy reordered eigenvalues to W.
    for (k = 0; k < n; k++) {
	w[k] = t[k + k * ldt];
    }
    work[1] = lwmin;
    return;
}
