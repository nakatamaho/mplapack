/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtrsna.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rtrsna(const char *job, const char *howmny, LOGICAL * select,
	    INTEGER n, REAL * t, INTEGER ldt, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, REAL * s, REAL * sep,
	    INTEGER mm, INTEGER * m, REAL * work, INTEGER ldwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, k, n2;
    REAL cs;
    INTEGER nn, ks;
    REAL sn, mu = 0.0, eps, est;
    INTEGER kase;
    REAL cond;
    INTEGER pair;
    INTEGER ierr;
    REAL dumm = 0.0, prod;
    INTEGER ifst;
    REAL lnrm;
    INTEGER ilst;
    REAL rnrm;
    REAL prod1, prod2, scale, delta;
    INTEGER isave[3];
    INTEGER wants;
    REAL dummy[1];
    REAL bignum;
    INTEGER wantbh;
    INTEGER somcon;
    REAL smlnum;
    INTEGER wantsp;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

//Decode and test the input parameters
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantsp = Mlsame(job, "V") || wantbh;
    somcon = Mlsame(howmny, "S");
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
    } else {
//Set M to the number of eigenpairs for which condition numbers
//are required, and test MM.
	if (somcon) {
	    *m = 0;
	    pair = MFALSE;
	    for (k = 0; k < n; k++) {
		if (pair) {
		    pair = MFALSE;
		} else {
		    if (k < n) {
			if (t[k + 1 + k * ldt] == Zero) {
			    if (select[k]) {
				++(*m);
			    }
			} else {
			    pair = MTRUE;
			    if (select[k] || select[k + 1]) {
				(*m) = (*m) + 2;
			    }
			}
		    } else {
			if (select[n]) {
			    ++(*m);
			}
		    }
		}
	    }
	} else {
	    (*m) = n;
	}
	if (mm < (*m)) {
	    *info = -13;
	} else if (ldwork < 1 || (wantsp && ldwork < n)) {
	    *info = -16;
	}
    }
    if (*info != 0) {
	Mxerbla("Rtrsna", -(*info));
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
    //Rlabad(&smlnum, &bignum);
    ks = 0;
    pair = MFALSE;
    for (k = 0; k < n; k++) {
//Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
	if (pair) {
	    pair = MFALSE;
	    goto L60;
	} else {
	    if (k < n) {
		pair = t[k + 1 + k * ldt] != Zero;
	    }
	}
//Determine whether condition numbers are required for the k-th
//eigenpair.
	if (somcon) {
	    if (pair) {
		if (!select[k] && !select[k + 1]) {
		    goto L60;
		}
	    } else {
		if (!select[k]) {
		    goto L60;
		}
	    }
	}
	ks++;
	if (wants) {
//Compute the reciprocal condition number of the k-th
//eigenvalue.
	    if (!pair) {
//Real eigenvalue.
		prod = Rdot(n, &vr[ks * ldvr + 1], 1, &vl[ks * ldvl + 1], 1);
		rnrm = Rnrm2(n, &vr[ks * ldvr + 1], 1);
		lnrm = Rnrm2(n, &vl[ks * ldvl + 1], 1);
		s[ks] = abs(prod) / (rnrm * lnrm);
	    } else {
//Complex eigenvalue.
		prod1 = Rdot(n, &vr[ks * ldvr + 1], 1, &vl[ks * ldvl + 1], 1);
		prod1 = prod1 + Rdot(n, &vr[(ks + 1) * ldvr + 1], 1, &vl[(ks + 1) * ldvl + 1], 1);
		prod2 = Rdot(n, &vl[ks * ldvl + 1], 1, &vr[(ks + 1) * ldvr + 1], 1);
		prod2 = prod2 - Rdot(n, &vl[(ks + 1) * ldvl + 1], 1, &vr[ks * ldvr + 1], 1);
		mtemp1 = Rnrm2(n, &vr[ks * ldvr + 1], 1);
		mtemp2 = Rnrm2(n, &vr[(ks + 1) * ldvr + 1], 1);
		rnrm = Rlapy2(mtemp1, mtemp2);
		mtemp1 = Rnrm2(n, &vl[ks * ldvl + 1], 1);
		mtemp2 = Rnrm2(n, &vl[(ks + 1) * ldvl + 1], 1);
		lnrm = Rlapy2(mtemp1, mtemp2);
		cond = Rlapy2(prod1, prod2) / (rnrm * lnrm);
		s[ks] = cond;
		s[ks + 1] = cond;
	    }
	}
	if (wantsp) {
//Estimate the reciprocal condition number of the k-th
//eigenvector.
//Copy the matrix T to the array WORK and swap the diagonal
//block beginning at T(k,k) to the (1,1) position.
	    Rlacpy("Full", n, n, &t[0], ldt, &work[0], ldwork);
	    ifst = k;
	    ilst = 1;
	    Rtrexc("No Q", n, &work[0], ldwork, dummy, 1, &ifst, &ilst, &work[(n + 1) * ldwork + 1], &ierr);
	    if (ierr == 1 || ierr == 2) {
//Could not swap because blocks not well separated
		scale = One;
		est = bignum;
	    } else {
//Reordering successful
		if (work[ldwork + 2] == Zero) {
//Form C = T22 - lambda*I in WORK(2:N,2:N).
		    for (i = 1; i < n; i++) {
			work[i + i * ldwork] = work[i + i * ldwork] - work[ldwork + 1];
		    }
		    n2 = 1;
		    nn = n - 1;
		} else {
/*                 Triangularize the 2 by 2 block by unitary */
/*                 transformation U = [  cs   i*ss ] */
/*                                    [ i*ss   cs  ]. */
/*                 such that the (1,1) position of WORK is complex */
/*                 eigenvalue lambda with positive imaginary part. (2,2) */
/*                 position of WORK is the complex eigenvalue lambda */
/*                 with negative imaginary  part. */
		    mu = sqrt(abs(work[(ldwork << 1) + 1])) * sqrt(abs(work[ldwork + 2]));
		    delta = Rlapy2(mu, work[ldwork + 2]);
		    cs = mu / delta;
		    sn = -work[ldwork + 2] / delta;
/*                 Form */
/*                 C' = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ] */
/*                                        [   mu                     ] */
/*                                        [         ..               ] */
/*                                        [             ..           ] */
/*                                        [                  mu      ] */
/*                 where C' is conjugate transpose of complex matrix C, */
/*                 and RWORK is stored starting in the N+1-st column of */
/*                 WORK. */
		    for (j = 3; j <= n; j++) {
			work[j * ldwork + 2] = cs * work[j * ldwork + 2];
			work[j + j * ldwork] = work[j + j * ldwork] - work[ldwork + 1];
		    }
		    work[(ldwork << 1) + 2] = Zero;
		    work[(n + 1) * ldwork + 1] = mu * Two;
		    for (i = 1; i < n - 1; i++) {
			work[i + (n + 1) * ldwork] = sn * work[(i + 1) * ldwork + 1];
		    }
		    n2 = 2;
		    nn = (n - 1) * 2;
		}
//Estimate norm(inv(C'))
		est = Zero;
		kase = 0;
	      L50:
		Rlacn2(nn, &work[(n + 2) * ldwork + 1], &work[(n + 4) * ldwork + 1], &iwork[1], &est, &kase, isave);
		if (kase != 0) {
		    if (kase == 1) {
			if (n2 == 1) {
//Real eigenvalue: solve C'*x = scale*c.
			    Rlaqtr(MTRUE, MTRUE, n - 1, &work[(ldwork << 1) + 2], ldwork, dummy, dumm, &scale, &work[(n + 4) * ldwork + 1], &work[(n + 6) * ldwork + 1], &ierr);
			} else {
//Complex eigenvalue: solve
//C'*(p+iq) = scale*(c+id) in real arithmetic.
			    Rlaqtr(MTRUE, MFALSE, n - 1, &work[(ldwork << 1) + 2], ldwork, &work[(n +
												  1) * ldwork + 1], mu, &scale,
				   &work[(n + 4) * ldwork + 1], &work[(n + 6) * ldwork + 1], &ierr);
			}
		    } else {
			if (n2 == 1) {
//Real eigenvalue: solve C*x = scale*c.
			    Rlaqtr(MFALSE, MTRUE, n - 1, &work[(ldwork << 1) + 2], ldwork, dummy, dumm, &scale, &work[(n + 4) * ldwork + 1], &work[(n + 6) * ldwork + 1], &ierr);
			} else {
//Complex eigenvalue: solve
//C*(p+iq) = scale*(c+id) in real arithmetic.
			    Rlaqtr(MFALSE, MFALSE, n - 1, &work[(ldwork << 1) + 2], ldwork, &work[(n +
												   1) * ldwork + 1], mu, &scale,
				   &work[(n + 4) * ldwork + 1], &work[(n + 6) * ldwork + 1], &ierr);

			}
		    }
		    goto L50;
		}
	    }
	    sep[ks] = scale / max(est, smlnum);
	    if (pair) {
		sep[ks + 1] = sep[ks];
	    }
	}
	if (pair) {
	    ks++;
	}
      L60:
	;
    }
    return;
}
