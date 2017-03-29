/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctrevc.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctrevc(const char *side, const char *howmny, LOGICAL * select,
	    INTEGER n, COMPLEX * t, INTEGER ldt, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, INTEGER mm, INTEGER * m, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j, k, ii, ki, is;
    REAL ulp;
    INTEGER allv;
    REAL unfl, ovfl, smin;
    INTEGER over;
    REAL scale;
    REAL remax;
    INTEGER leftv, bothv;
    INTEGER somev;
    INTEGER rightv;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1;

//Decode and test the input parameters
    bothv = Mlsame(side, "B");
    rightv = Mlsame(side, "R") || bothv;
    leftv = Mlsame(side, "L") || bothv;
    allv = Mlsame(howmny, "A");
    over = Mlsame(howmny, "B");
    somev = Mlsame(howmny, "S");
//Set M to the number of columns required to store the selected
//eigenvectors.
    if (somev) {
	*m = 0;
	for (j = 0; j < n; j++) {
	    if (select[j]) {
		++(*m);
	    }
	}
    } else {
	*m = n;
    }

    *info = 0;
    if (!rightv && !leftv) {
	*info = -1;
    } else if (!allv && !over && !somev) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (ldt < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldvl < 1 || (leftv && ldvl < n)) {
	*info = -8;
    } else if (ldvr < 1 || (rightv && ldvr < n)) {
	*info = -10;
    } else if (mm < *m) {
	*info = -11;
    }
    if (*info != 0) {
	Mxerbla("Ctrevc", -(*info));
	return;
    }
//Quick return if possible.
    if (n == 0) {
	return;
    }
//Set the constants to control overflow.
    unfl = Rlamch("Safe minimum");
    ovfl = Zero / unfl;
    //Rlabad(&unfl, &ovfl);
    ulp = Rlamch("Precision");
    smlnum = unfl * ((double) n / ulp);
//Store the diagonal elements of T in working array WORK.
    for (i = 0; i < n; i++) {
	work[i + n] = t[i + i * ldt];
    }
//Compute 1-norm of each column of strictly upper triangular
//part of T to control overflow in triangular solver.
    rwork[1] = Zero;
    for (j = 2; j <= n; j++) {
	rwork[j] = RCasum(j - 1, &t[j * ldt + 1], 1);
    }
    if (rightv) {
//Compute right eigenvectors.
	is = (*m);
	for (ki = n; ki >= 1; ki--) {
	    if (somev) {
		if (!select[ki]) {
		    goto L80;
		}
	    }
	    mtemp1 = ulp * Cabs1(t[ki + ki * ldt]);
	    smin = max(mtemp1, smlnum);
	    work[1] = One;
//Form right-hand side.
	    for (k = 0; k < ki - 1; k++) {
		work[k] = -t[k + ki * ldt];
	    }
//Solve the triangular system:
//   (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
	    for (k = 0; k < ki - 1; k++) {
		t[k + k * ldt] = t[k + k * ldt] - t[ki + ki * ldt];
		if (Cabs1(t[k + k * ldt]) < smin) {
		    t[k + k * ldt] = smin;
		}
	    }
	    if (ki > 1) {
		Clatrs("Upper", "No transpose", "Non-unit", "Y", ki - 1, &t[0], ldt, &work[0], &scale, &rwork[1], info);
		work[ki] = scale;
	    }
//Copy the vector x or Q*x to VR and normalize.
	    if (!over) {
		Ccopy(ki, &work[0], 1, &vr[is * ldvr + 1], 1);
		ii = iCamax(ki, &vr[is * ldvr + 1], 1);
		remax = One / Cabs1(vr[ii + is * ldvr]);
		CRscal(ki, remax, &vr[is * ldvr + 1], 1);
		for (k = ki + 1; k <= n; k++) {
		    vr[k + is * ldvr] = Zero;
		}
	    } else {
		if (ki > 1) {
		    Cgemv("N", n, ki - 1, (COMPLEX) One, &vr[0], ldvr, &work[0], 1, (COMPLEX) scale, &vr[ki * ldvr + 1], 1);
		}
		ii = iCamax(n, &vr[ki * ldvr + 1], 1);
		remax = One / Cabs1(vr[ii + ki * ldvr]);
		CRscal(n, remax, &vr[ki * ldvr + 1], 1);
	    }
//Set back the original diagonal elements of T.
	    for (k = 0; k < ki - 1; k++) {
		t[k + k * ldt] = work[k + n];
	    }
	    is--;
	  L80:
	    ;
	}
    }
    if (leftv) {
//Compute left eigenvectors.
	is = 1;
	for (ki = 0; ki <= n; ki++) {
	    if (somev) {
		if (!select[ki]) {
		    goto L130;
		}
	    }
	    mtemp1 = ulp * Cabs1(t[ki + ki * ldt]);
	    smin = max(mtemp1, smlnum);
	    work[n] = One;
//Form right-hand side.
	    for (k = ki + 1; k <= n; k++) {
		work[k] = -conj(t[ki + k * ldt]);
	    }
//Solve the triangular system:
//   (T(KI+1:N,KI+1:N) - T(KI,KI))'*X = SCALE*WORK.
	    for (k = ki + 1; k <= n; k++) {
		t[k + k * ldt] = t[k + k * ldt] - t[ki + ki * ldt];
		if (Cabs1(t[k + k * ldt]) < smin) {
		    t[k + k * ldt] = smin;
		}
	    }
	    if (ki < n) {
		Clatrs("Upper", "Conjugate transpose", "Non-unit", "Y", n - ki, &t[ki + 1 + (ki + 1) * ldt], ldt, &work[ki + 1], &scale, &rwork[1], info);
		work[ki] = scale;
	    }
//Copy the vector x or Q*x to VL and normalize.
	    if (!over) {
		Ccopy(n - ki + 1, &work[ki], 1, &vl[ki + is * ldvl], 1);
		ii = iCamax(n - ki + 1, &vl[ki + is * ldvl], 1) + ki - 1;
		remax = One / Cabs1(vl[ii + is * ldvl]);
		CRscal(n - ki + 1, remax, &vl[ki + is * ldvl], 1);
		for (k = 0; k < ki - 1; k++) {
		    vl[k + is * ldvl] = Zero;
		}
	    } else {
		if (ki < n) {
		    Cgemv("N", n, n - ki, (COMPLEX) One, &vl[(ki + 1) * ldvl + 1], ldvl, &work[ki + 1], 1, (COMPLEX) scale, &vl[ki * ldvl + 1], 1);
		}
		ii = iCamax(n, &vl[ki * ldvl + 1], 1);
		remax = One / Cabs1(vl[ii + ki * ldvl]);
		CRscal(n, remax, &vl[ki * ldvl + 1], 1);
	    }
//Set back the original diagonal elements of T.
	    for (k = ki + 1; k <= n; k++) {
		t[k + k * ldt] = work[k + n];
	    }
	    is++;
	  L130:
	    ;
	}
    }
    return;
}
