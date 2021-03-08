/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtrsen.cpp,v 1.12 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void
Rtrsen(const char *job, const char *compq, LOGICAL * select, INTEGER
       n, REAL * t, INTEGER ldt, REAL * q, INTEGER ldq,
       REAL * wr, REAL * wi, INTEGER m, REAL * s, REAL * sep, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER k, n1, n2, kk, nn, ks;
    REAL est;
    INTEGER kase;
    INTEGER pair;
    INTEGER ierr;
    INTEGER swap;
    REAL scale;
    INTEGER isave[3], lwmin = 0;
    INTEGER wantq, wants;
    REAL rnorm;
    INTEGER wantbh;
    INTEGER liwmin = 0;
    INTEGER wantsp, lquery;
    REAL Zero = 0.0, One = 1.0;

//Decode and test the input parameters
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantsp = Mlsame(job, "V") || wantbh;
    wantq = Mlsame(compq, "V");

    *info = 0;
    lquery = lwork == -1;
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
    } else {
//Set M to the dimension of the specified invariant subspace,
//and test LWORK and LIWORK.
	m = 0;
	pair = MFALSE;
	for (k = 0; k < n; k++) {
	    if (pair) {
		pair = MFALSE;
	    } else {
		if (k < n) {
		    if (t[k + 1 + k * ldt] == Zero) {
			if (select[k]) {
			    ++(m);
			}
		    } else {
			pair = MTRUE;
			if (select[k] || select[k + 1]) {
			    m = m + 2;
			}
		    }
		} else {
		    if (select[n]) {
			++(m);
		    }
		}
	    }
	}
	n1 = m;
	n2 = n - m;
	nn = n1 * n2;
	if (wantsp) {
	    lwmin = max((INTEGER) 1, nn * 2);
	    liwmin = max((INTEGER) 1, nn);
	} else if (Mlsame(job, "N")) {
	    lwmin = max((INTEGER) 1, n);
	    liwmin = 1;
	} else if (Mlsame(job, "E")) {
	    lwmin = max((INTEGER) 1, nn);
	    liwmin = 1;
	}
	if (lwork < lwmin && !lquery) {
	    *info = -15;
	} else if (liwork < liwmin && !lquery) {
	    *info = -17;
	}
    }
    if (*info == 0) {
	work[1] = lwmin;
	iwork[1] = liwmin;
    }
    if (*info != 0) {
	Mxerbla("Rtrsen", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible.
    if (m == n || m == 0) {
	if (wants) {
	    *s = One;
	}
	if (wantsp) {
	    *sep = Rlange("1", n, n, t, ldt, work);
	}
	goto L40;
    }
//Collect the selected blocks at the top-left corner of T.
    ks = 0;
    pair = MFALSE;
    for (k = 0; k < n; k++) {
	if (pair) {
	    pair = MFALSE;
	} else {
	    swap = select[k];
	    if (k < n) {
		if (t[k + 1 + k * ldt] != Zero) {
		    pair = MTRUE;
		    swap = swap || select[k + 1];
		}
	    }
	    if (swap) {
		ks++;
//Swap the K-th block to position KS.
		ierr = 0;
		kk = k;
		if (k != ks) {
		    Rtrexc(compq, n, t, ldt, q, ldq, &kk, &ks, work, &ierr);
		}
		if (ierr == 1 || ierr == 2) {
//Blocks too close to swap: exit.
		    *info = 1;
		    if (wants) {
			*s = Zero;
		    }
		    if (wantsp) {
			*sep = Zero;
		    }
		    goto L40;
		}
		if (pair) {
		    ks++;
		}
	    }
	}
    }
    if (wants) {
//Solve Sylvester equation for R:
//   T11*R - R*T22 = scale*T12
	Rlacpy("F", n1, n2, &t[(n1 + 1) * ldt + 1], ldt, work, n1);
	Rtrsyl("N", "N", -1, n1, n2, t, ldt, &t[n1 + 1 + (n1 + 1) * ldt], ldt, work, n1, &scale, &ierr);
//Estimate the reciprocal of the condition number of the cluster
//of eigenvalues.
	rnorm = Rlange("F", n1, n2, work, n1, work);
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
	Rlacn2(nn, &work[nn + 1], work, iwork, &est, &kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
//Solve  T11*R - R*T22 = scale*X.
		Rtrsyl("N", "N", -1, n1, n2, t, ldt, &t[n1 + 1 + (n1 + 1) * ldt], ldt, work, n1, &scale, &ierr);
	    } else {
//Solve  T11'*R - R*T22' = scale*X.
		Rtrsyl("T", "T", -1, n1, n2, t, ldt, &t[n1 + 1 + (n1 + 1) * ldt], ldt, work, n1, &scale, &ierr);
	    }
	    goto L30;
	}
	*sep = scale / est;
    }
  L40:
//store the output eigenvalues in wr and wi.
    for (k = 0; k < n; k++) {
	wr[k] = t[k + k * ldt];
	wi[k] = Zero;
    }
    for (k = 0; k < n - 1; k++) {
	if (t[k + 1 + k * ldt] != Zero) {
	    wi[k] = sqrt(abs(t[k + (k + 1) * ldt])) * sqrt(abs(t[k + 1 + k * ldt]));
	    wi[k + 1] = -wi[k];
	}
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
