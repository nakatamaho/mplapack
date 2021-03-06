/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chsein.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chsein(const char *side, const char *eigsrc, const char *initv, LOGICAL *
	    select, INTEGER n, COMPLEX * h, INTEGER ldh, COMPLEX *
	    w, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, INTEGER mm, INTEGER * m, COMPLEX * work, REAL * rwork, INTEGER * ifaill, INTEGER * ifailr, INTEGER * info)
{
    INTEGER i, k, kl, kr, ks;
    COMPLEX wk;
    INTEGER kln;
    REAL ulp, eps3 = 0, unfl;
    INTEGER iinfo;
    INTEGER leftv, bothv;
    REAL hnorm;
    INTEGER noinit;
    INTEGER ldwork;
    INTEGER rightv, fromqr;
    REAL smlnum;
    REAL Zero = 0.0;

//Decode and test the input parameters.
    bothv = Mlsame(side, "B");
    rightv = Mlsame(side, "R") || bothv;
    leftv = Mlsame(side, "L") || bothv;
    fromqr = Mlsame(eigsrc, "Q");
    noinit = Mlsame(initv, "N");
//Set M to the number of columns required to store the selected
//eigenvectors.
    (*m) = 0;
    for (k = 0; k < n; k++) {
	if (select[k]) {
	    ++(*m);
	}
    }
    *info = 0;
    if (!rightv && !leftv) {
	*info = -1;
    } else if (!fromqr && !Mlsame(eigsrc, "N")) {
	*info = -2;
    } else if (!noinit && !Mlsame(initv, "U")) {
	*info = -3;
    } else if (n < 0) {
	*info = -5;
    } else if (ldh < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldvl < 1 || (leftv && ldvl < n)) {
	*info = -10;
    } else if (ldvr < 1 || (rightv && ldvr < n)) {
	*info = -12;
    } else if (mm < (*m)) {
	*info = -13;
    }
    if (*info != 0) {
	Mxerbla("Chsein", -(*info));
	return;
    }
//Quick return if possible.
    if (n == 0) {
	return;
    }
//Set machine-dependent constants.
    unfl = Rlamch("Safe minimum");
    ulp = Rlamch("Precision");
    smlnum = unfl * (n / ulp);
    ldwork = n;
    kl = 0;
    kln = 0;
    if (fromqr) {
	kr = 0;
    } else {
	kr = n;
    }
    ks = 1;
    for (k = 0; k < n; k++) {
	if (select[k]) {
//Compute eigenvector(s) corresponding to W(K).
	    if (fromqr) {
//If affiliation of eigenvalues is known, check whether
//the matrix splits.
//Determine KL and KR such that 1 <= KL <= K <= KR <= N
//and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
//KR = N).
//Then inverse iteration can be performed with the
//submatrix H(KL:N,KL:N) for a left eigenvector, and with
//the submatrix H(1:KR,1:KR) for a right eigenvector.
		for (i = k; i >= kl + 1; i--) {
		    if (h[i + (i - 1) * ldh] == Zero) {
			goto L30;
		    }
		}
	      L30:
		kl = i;
		if (k > kr) {
		    for (i = k; i <= n - 1; i++) {
			if (h[i + 1 + i * ldh] == Zero) {
			    goto L50;
			}
		    }
		  L50:
		    kr = i;
		}
	    }
	    if (kl != kln) {
		kln = kl;
//Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
//has not ben computed before.
		hnorm = Clanhs("I", kr - kl + 1, &h[kl + kl * ldh], ldh, &rwork[1]);
		if (hnorm > Zero) {
		    eps3 = hnorm * ulp;
		} else {
		    eps3 = smlnum;
		}
	    }
//Perturb eigenvalue if it is close to any previous
//selected eigenvalues affiliated to the submatrix
//H(KL:KR,KL:KR). Close roots are modified by EPS3.
	    wk = w[k];
	  L60:
	    for (i = k - 1; i >= kl; i--) {
		if (select[i] && (abs((w[i] - wk).real()) + abs((w[i] - wk).imag())) < eps3) {
		    wk = wk + eps3;
		    goto L60;
		}
	    }
	    w[k] = wk;
	    if (leftv) {
//Compute left eigenvector.
		Claein(MFALSE, noinit, n - kl + 1, &h[kl + kl * ldh], ldh, wk, &vl[kl + ks * ldvl], &work[0], ldwork, &rwork[1], eps3, smlnum, &iinfo);
		if (iinfo > 0) {
		    ++(*info);
		    ifaill[ks] = k;
		} else {
		    ifaill[ks] = 0;
		}
		for (i = 0; i < kl - 1; i++) {
		    vl[i + ks * ldvl] = Zero;
		}
	    }
	    if (rightv) {
//Compute right eigenvector.
		Claein(MTRUE, noinit, kr, &h[0], ldh, wk, &vr[ks * ldvr + 1], &work[0], ldwork, &rwork[1], eps3, smlnum, &iinfo);
		if (iinfo > 0) {
		    ++(*info);
		    ifailr[ks] = k;
		} else {
		    ifailr[ks] = 0;
		}
		for (i = kr + 1; i <= n; i++) {
		    vr[i + ks * ldvr] = Zero;
		}
	    }
	    ks++;
	}
    }
    return;
}
