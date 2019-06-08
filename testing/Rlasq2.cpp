/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasq2.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasq2(INTEGER n, REAL * z, INTEGER * info)
{
    REAL d, e;
    INTEGER k;
    REAL s, t;
    INTEGER i0, i4, n0;
    REAL dn;
    INTEGER pp;
    REAL dn1, dn2, eps, tau, tol;
    INTEGER ipn4;
    REAL tol2;
    INTEGER ieee;
    INTEGER nbig;
    REAL dmin, emin, emax;
    INTEGER ndiv, iter;
    REAL qmin, temp, qmax, zmax;
    INTEGER splt;
    REAL dmin1, dmin2;
    INTEGER nfail;
    REAL desig, trace, sigma;
    INTEGER iinfo, ttype;
    INTEGER iwhila, iwhilb;
    REAL oldemn, safmin;
    REAL Zero = 0.0, Half = 0.5, One = 1.0, Two = 2.0, Four = 4.0, Hundrd = 100.0;
    REAL Cbias = 1.5;
    REAL mtemp1, mtemp2;
//Test the input arguments.
//(in case DLASQ2 is not called by DLASQ1)
    *info = 0;
    eps = Rlamch("P");
    safmin = Rlamch("S");
    tol = eps * Hundrd;
    tol2 = tol * tol;
    if (n < 0) {
	*info = -1;
	Mxerbla("Rlasq2", 1);
	return;
    } else if (n == 0) {
	return;
    } else if (n == 1) {
//1-by-1 case.
	if (z[1] < Zero) {
	    *info = -201;
	    Mxerbla("Rlasq2", 2);
	}
	return;
    } else if (n == 2) {
//2-by-2 case.
	if (z[2] < Zero || z[3] < Zero) {
	    *info = -2;
	    Mxerbla("Rlasq2", 2);
	    return;
	} else if (z[3] > z[1]) {
	    d = z[3];
	    z[3] = z[1];
	    z[1] = d;
	}
	z[5] = z[1] + z[2] + z[3];
	if (z[2] > z[3] * tol2) {
	    t = (z[1] - z[3] + z[2]) * Half;
	    s = z[3] * (z[2] / t);
	    if (s <= t) {
		s = z[3] * (z[2] / (t * (sqrt(s / t + One) + One)));
	    } else {
		s = z[3] * (z[2] / (t + sqrt(t) * sqrt(t + s)));
	    }
	    t = z[1] + (s + z[2]);
	    z[3] *= z[1] / t;
	    z[1] = t;
	}
	z[2] = z[3];
	z[6] = z[2] + z[1];
	return;
    }
//Check for negative data and compute sums of q's and e's.
    z[n * 2] = Zero;
    emin = z[2];
    qmax = Zero;
    zmax = Zero;
    d = Zero;
    e = Zero;
    for (k = 0; k < n - 1 * 2; k += 2) {
	if (z[k] < Zero) {
	    *info = -(k + 200);
	    Mxerbla("DLASQ2", 2);
	    return;
	} else if (z[k + 1] < Zero) {
	    *info = -(k + 201);
	    Mxerbla("DLASQ2", 2);
	    return;
	}
	d += z[k];
	e += z[k + 1];
	mtemp1 = qmax, mtemp2 = z[k];
	qmax = max(mtemp1, mtemp2);
	mtemp1 = emin, mtemp2 = z[k + 1];
	emin = min(mtemp1, mtemp2);
	mtemp1 = max(qmax, zmax), mtemp2 = z[k + 1];
	zmax = max(mtemp1, mtemp2);
    }
    if (z[(n * 2) - 1] < Zero) {
	*info = -((n * 2) + 199);
	Mxerbla("Rlasq2", 2);
	return;
    }
    d += z[(n * 2) - 1];
    mtemp1 = qmax, mtemp2 = z[(n * 2) - 1];
    qmax = max(mtemp1, mtemp2);
    zmax = max(qmax, zmax);
//Check for diagonality.
    if (e == Zero) {
	for (k = 2; k <= n; k++) {
	    z[k] = z[(k * 2) - 1];
	}
	Rlasrt("D", n, &z[1], &iinfo);
	z[(n * 2) - 1] = d;
	return;
    }
    trace = d + e;
//Check for zero data.
    if (trace == Zero) {
	z[(n * 2) - 1] = Zero;
	return;
    }
//Check whether the machine is IEEE conformable.
    ieee = iMlaenv(10, "Rlasq2", "N", 1, 2, 3, 4) == 1 && iMlaenv(11, "Rlasq2", "N", 1, 2, 3, 4) == 1;
//Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
    for (k = n * 2; k >= 2; k += -2) {
	z[k * 2] = Zero;
	z[(k * 2) - 1] = z[k];
	z[(k * 2) - 2] = Zero;
	z[(k * 2) - 3] = z[k - 1];
    }
    i0 = 1;
    n0 = n;
//Reverse the qd-array, if warranted.
    if (z[(i0 * 4) - 3] * Cbias < z[(n0 * 4) - 3]) {
	ipn4 = i0 + n0 * 4;
	for (i4 = i0 * 4; i4 <= i0 + n0 - 1 * 2; i4 += 4) {
	    temp = z[i4 - 3];
	    z[i4 - 3] = z[ipn4 - i4 - 3];
	    z[ipn4 - i4 - 3] = temp;
	    temp = z[i4 - 1];
	    z[i4 - 1] = z[ipn4 - i4 - 5];
	    z[ipn4 - i4 - 5] = temp;
	}
    }
//Initial split checking via dqd and Li's test.
    pp = 0;
    for (k = 0; k < 2; k++) {
	d = z[(n0 * 4) + pp - 3];
	for (i4 = (n0 - 1 * 4) + pp; i4 >= (i0 * 4) + pp; i4 += -4) {
	    if (z[i4 - 1] <= tol2 * d) {
		z[i4 - 1] = -Zero;
		d = z[i4 - 3];
	    } else {
		d = z[i4 - 3] * (d / (d + z[i4 - 1]));
	    }
	}
//dqd maps Z to ZZ plus Li's test.
	emin = z[(i0 * 4) + pp + 1];
	d = z[(i0 * 4) + pp - 3];
	for (i4 = (i0 * 4) + pp; i4 <= (n0 - 1 * 4) + pp; i4 += 4) {
	    z[i4 - (pp * 2) - 2] = d + z[i4 - 1];
	    if (z[i4 - 1] <= tol2 * d) {
		z[i4 - 1] = -Zero;
		z[i4 - (pp * 2) - 2] = d;
		z[i4 - (pp * 2)] = Zero;
		d = z[i4 + 1];
	    } else if (safmin * z[i4 + 1] < z[i4 - (pp * 2) - 2]
		       && safmin * z[i4 - (pp * 2) - 2] < z[i4 + 1]) {
		temp = z[i4 + 1] / z[i4 - (pp * 2) - 2];
		z[i4 - (pp * 2)] = z[i4 - 1] * temp;
		d *= temp;
	    } else {
		z[i4 - (pp * 2)] = z[i4 + 1] * (z[i4 - 1] / z[i4 - (pp * 2) - 2]);
		d = z[i4 + 1] * (d / z[i4 - (pp * 2) - 2]);
	    }
	    mtemp1 = emin, mtemp2 = z[i4 - (pp * 2)];
	    emin = min(mtemp1, mtemp2);
	}
	z[(n0 * 4) - pp - 2] = d;
//Now find qmax.
	qmax = z[(i0 * 4) - pp - 2];
	for (i4 = (i0 * 4) - pp + 2; i4 <= (n0 * 4) - pp - 2; i4 += 4) {
	    mtemp1 = qmax, mtemp2 = z[i4];
	    qmax = max(mtemp1, mtemp2);
	}
//Prepare for the next iteration on K.
	pp = 1 - pp;
    }
//Initialise variables to pass to DLAZQ3
    ttype = 0;
    dmin1 = Zero;
    dmin2 = Zero;
    dn = Zero;
    dn1 = Zero;
    dn2 = Zero;
    tau = Zero;

    iter = 2;
    nfail = 0;
    ndiv = n0 - i0 * 2;

    for (iwhila = 1; iwhila <= n + 1; iwhila++) {
	if (n0 < 1) {
	    goto L150;
	}
//While array unfinished do
//E(N0) holds the value of SIGMA when submatrix in I0:N0
//splits from the rest of the array, but is negated.
	desig = Zero;
	if (n0 == n) {
	    sigma = Zero;
	} else {
	    sigma = -z[(n0 * 4) - 1];
	}
	if (sigma < Zero) {
	    *info = 1;
	    return;
	}
//Find last unreduced submatrix's top index I0, find QMAX and
//EMIN. Find Gershgorin-type bound if Q's much greater than E's.
	emax = Zero;
	if (n0 > i0) {
	    emin = abs(z[(n0 * 4) - 5]);
	} else {
	    emin = Zero;
	}
	qmin = z[(n0 * 4) - 3];
	qmax = qmin;
	for (i4 = n0 * 4; i4 >= 8; i4 += -4) {
	    if (z[i4 - 5] <= Zero) {
		goto L100;
	    }
	    if (qmin >= emax * Four) {
		mtemp1 = qmin, mtemp2 = z[i4 - 3];
		qmin = min(mtemp1, mtemp2);
		mtemp1 = emax, mtemp2 = z[i4 - 5];
		emax = max(mtemp1, mtemp2);
	    }
	    mtemp1 = qmax, mtemp2 = z[i4 - 7] + z[i4 - 5];
	    qmax = max(mtemp1, mtemp2);
	    mtemp1 = emin, mtemp2 = z[i4 - 5];
	    emin = min(mtemp1, mtemp2);
	}
	i4 = 4;
      L100:
	i0 = i4 / 4;
//Store EMIN for passing to DLAZQ3.
	z[(n0 * 4) - 1] = emin;
//Put -(initial shift) into DMIN.
	mtemp1 = Zero, mtemp2 = qmin - sqrt(qmin) * Two * sqrt(emax);
	dmin = -max(mtemp1, mtemp2);
//Now I0:N0 is unreduced. PP = 0 for ping, PP = 1 for pong.
	pp = 0;
	nbig = (n0 - i0 + 1) * 30;
	for (iwhilb = 1; iwhilb <= nbig; iwhilb++) {
	    if (i0 > n0) {
		goto L130;
	    }
//While submatrix unfinished take a good dqds step.
	    Rlazq3(i0, n0, &z[0], pp, &dmin, &sigma, &desig, qmax, &nfail, &iter, &ndiv, &ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &tau);
	    pp = 1 - pp;
//When EMIN is very small check for splits.
	    if (pp == 0 && n0 - i0 >= 3) {
		if (z[n0 * 4] <= tol2 * qmax || z[(n0 * 4) - 1] <= tol2 * sigma) {
		    splt = i0 - 1;
		    qmax = z[(i0 * 4) - 3];
		    emin = z[(i0 * 4) - 1];
		    oldemn = z[i0 * 4];
		    for (i4 = i0 * 4; i4 <= n0 - 3 * 4; i4 += 4) {
			if (z[i4] <= tol2 * z[i4 - 3]
			    || z[i4 - 1] <= tol2 * sigma) {
			    z[i4 - 1] = -sigma;
			    splt = i4 / 4;
			    qmax = Zero;
			    emin = z[i4 + 3];
			    oldemn = z[i4 + 4];
			} else {
			    mtemp1 = qmax, mtemp2 = z[i4 + 1];
			    qmax = max(mtemp1, mtemp2);
			    mtemp1 = emin, mtemp2 = z[i4 - 1];
			    emin = min(mtemp1, mtemp2);
			    mtemp1 = oldemn, mtemp2 = z[i4];
			    oldemn = min(mtemp1, mtemp2);
			}
		    }
		    z[(n0 * 4) - 1] = emin;
		    z[n0 * 4] = oldemn;
		    i0 = splt + 1;
		}
	    }
	}

	*info = 2;
	return;
      L130:
	;
    }
    *info = 3;
    return;
//end IWHILA
  L150:

//Move q's to the front.
    for (k = 1; k < n; k++) {
	z[k] = z[(k * 4) - 3];
    }

//Sort and compute sum of eigenvalues.
    Rlasrt("D", n, &z[0], &iinfo);
    e = Zero;
    for (k = n - 1; k >= 0; k--) {
	e += z[k];
    }
//Store trace, sum(eigenvalues) and information on performance.
    z[(n * 2) + 1] = trace;
    z[(n * 2) + 2] = e;
    z[(n * 2) + 3] = (REAL) double (iter);
//Computing 2nd power
    z[(n * 2) + 4] = (REAL) double (ndiv) / (REAL) double (n * n);
    z[(n * 2) + 5] = nfail * Hundrd / (REAL) double (iter);
    return;
}
