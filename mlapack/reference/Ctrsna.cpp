/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctrsna.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctrsna(const char *job, const char *howmny, LOGICAL * select,
	    INTEGER n, COMPLEX * t, INTEGER ldt, COMPLEX * vl,
	    INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, REAL * s, REAL * sep, INTEGER mm, INTEGER * m, COMPLEX * work, INTEGER ldwork, REAL * rwork, INTEGER * info)
{
    INTEGER i, j, k, ks, ix;
    REAL eps, est;
    INTEGER kase, ierr;
    COMPLEX prod;
    REAL lnrm, rnrm, scale;
    INTEGER isave[3];
    COMPLEX dummy[1];
    INTEGER wants;
    REAL xnorm;
    REAL bignum;
    INTEGER wantbh;
    INTEGER somcon;
    char normin;
    REAL smlnum;
    INTEGER wantsp;
    REAL Zero = 0.0, One = 1.0;

//Decode and test the input parameters
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantsp = Mlsame(job, "V") || wantbh;
    somcon = Mlsame(howmny, "S");
//Set M to the number of eigenpairs for which condition numbers are
//to be computed.
    if (somcon) {
	(*m) = 0;
	for (j = 0; j < n; j++) {
	    if (select[j]) {
		++(*m);
	    }
	}
    } else {
	(*m) = n;
    }
    *info = 0;
    if (!wants && !wantsp) {
	*info = -1;
    } else if (!Mlsame(howmny, "A") && !somcon) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (ldt < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldvl < 1 || (wants && ldvl < n)) {
	*info = -8;
    } else if (ldvr < 1 || (wants && ldvr < n)) {
	*info = -10;
    } else if (mm < (*m)) {
	*info = -13;
    } else if (ldwork < 1 || (wantsp && ldwork < n)) {
	*info = -16;
    }
    if (*info != 0) {
	Mxerbla("Ctrnsa", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (somcon) {
	    if (!select[1]) {
		return;
	    }
	}
	if (wants) {
	    s[1] = One;
	}
	if (wantsp) {
	    sep[1] = abs(t[ldt + 1]);
	}
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = One / smlnum;
    ks = 1;
    for (k = 0; k < n; k++) {
	if (somcon) {
	    if (!select[k]) {
		goto L50;
	    }
	}
	if (wants) {
//Compute the reciprocal condition number of the k-th
//eigenvalue.
	    prod = Cdotc(n, &vr[ks * ldvr + 1], 1, &vl[ks * ldvl + 1], 1);
	    rnrm = RCnrm2(n, &vr[ks * ldvr + 1], 1);
	    lnrm = RCnrm2(n, &vl[ks * ldvl + 1], 1);
	    s[ks] = abs(prod) / (rnrm * lnrm);
	}
	if (wantsp) {
//Estimate the reciprocal condition number of the k-th
//eigenvector.
//Copy the matrix T to the array WORK and swap the k-th
//diagonal element to the (1,1) position.
	    Clacpy("Full", n, n, &t[0], ldt, &work[0], ldwork);
	    Ctrexc("No Q", n, &work[0], ldwork, &dummy[0], 1, k, 1, &ierr);
//Form  C = T22 - lambda*I in WORK(2:N,2:N).
	    for (i = 1; i < n; i++) {
		work[i + i * ldwork] = work[i + i * ldwork] - work[ldwork + 1];
	    }
//Estimate a lower bound for the 1-norm of inv(C'). The 1st
//and (N+1)th columns of WORK are used to store work vectors.
	    sep[ks] = Zero;
	    est = Zero;
	    kase = 0;
	    normin = 'N';
	  L30:
	    Clacn2(n - 1, &work[(n + 1) * ldwork + 1], &work[0], &est, &kase, isave);
	    if (kase != 0) {
		if (kase == 1) {
//Solve C'*x = scale*b
		    Clatrs("Upper", "Conjugate transpose", "Nonunit", (const char *) normin, n - 1, &work[(ldwork << 1) + 2], ldwork, &work[0], &scale, &rwork[1], &ierr);
		} else {
//Solve C*x = scale*b
		    Clatrs("Upper", "No transpose", "Nonunit", (const char *) normin, n - 1, &work[(ldwork << 1) + 2], ldwork, &work[0], &scale, &rwork[1], &ierr);
		}
		normin = 'Y';
		if (scale != One) {
//Multiply by 1/SCALE if doing so will not cause
//overflow.
		    ix = iCamax(n - 1, &work[0], 1);
		    xnorm = Cabs1(work[ix + ldwork]);
		    if (scale < xnorm * smlnum || scale == Zero) {
			goto L40;
		    }
		    CRrscl(n, scale, &work[0], 1);
		}
		goto L30;
	    }
	    sep[ks] = One / max(est, smlnum);
	}
      L40:
	ks++;
      L50:
	;
    }
    return;
}
