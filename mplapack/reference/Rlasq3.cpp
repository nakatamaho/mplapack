/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasq3.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasq3(INTEGER i0, INTEGER n0, REAL * z, INTEGER pp, REAL * dmin, REAL * sigma, REAL * desig, REAL qmax, INTEGER * nfail, INTEGER * iter, INTEGER * ndiv, INTEGER ieee)
{
    REAL s, t;
    INTEGER j4, nn;
    REAL eps, tol;
    INTEGER n0in, ipn4;
    REAL tol2, temp;
    REAL safmin;
    REAL Zero = 0.0, Qurtr = 0.25, Half = 0.5, One = 1.0, Two = 2.0, Hundrd = 100.0;
    REAL Cbias = 1.5, tau;

    REAL mtemp1, mtemp2, mtemp3;
    INTEGER ttype = 0;
    REAL dmin1 = Zero;
    REAL dmin2 = Zero;
    REAL dn = Zero;
    REAL dn1 = Zero;
    REAL dn2 = Zero;

    n0in = n0;
    eps = Rlamch("Precision");
    safmin = Rlamch("Safe minimum");
    tol = eps * Hundrd;
    tol2 = tol * tol;
//Check for deflation.
  L10:
    if (n0 < i0) {
	return;
    }
    if (n0 == i0) {
	goto L20;
    }
    nn = (n0 << 2) + pp;
    if (n0 == i0 + 1) {
	goto L40;
    }
//Check whether E(N0-1) is negligible, 1 eigenvalue.
    if (z[nn - 5] > tol2 * (*sigma + z[nn - 3])
	&& z[nn - (pp * 2) - 4] > tol2 * z[nn - 7]) {
	goto L30;
    }
  L20:
    z[(n0 << 2) - 3] = z[(n0 << 2) + pp - 3] + *sigma;
    --(n0);
    goto L10;
//Check  whether E(N0-2) is negligible, 2 eigenvalues.
  L30:
    if (z[nn - 9] > tol2 * *sigma && z[nn - (pp * 2) - 8] > tol2 * z[nn - 11]) {
	goto L50;
    }
  L40:
    if (z[nn - 3] > z[nn - 7]) {
	s = z[nn - 3];
	z[nn - 3] = z[nn - 7];
	z[nn - 7] = s;
    }
    if (z[nn - 5] > z[nn - 3] * tol2) {
	t = (z[nn - 7] - z[nn - 3] + z[nn - 5]) * Half;
	s = z[nn - 3] * (z[nn - 5] / t);
	if (s <= t) {
	    s = z[nn - 3] * (z[nn - 5] / (t * (sqrt(s / t + One) + One)));
	} else {
	    s = z[nn - 3] * (z[nn - 5] / (t + sqrt(t) * sqrt(t + s)));
	}
	t = z[nn - 7] + (s + z[nn - 5]);
	z[nn - 3] *= z[nn - 7] / t;
	z[nn - 7] = t;
    }
    z[(n0 << 2) - 7] = z[nn - 7] + *sigma;
    z[(n0 << 2) - 3] = z[nn - 3] + *sigma;
    n0 += -2;
    goto L10;
  L50:
//Reverse the qd-array, if warranted.
    if (*dmin <= Zero || n0 < n0in) {
	if (z[(i0 * 4) + pp - 3] * Cbias < z[(n0 * 4) + pp - 3]) {
	    ipn4 = i0 + n0 * 4;
	    for (j4 = i0 * 4; j4 <= i0 + n0 - 1 * 2; j4 += 4) {
		temp = z[j4 - 3];
		z[j4 - 3] = z[ipn4 - j4 - 3];
		z[ipn4 - j4 - 3] = temp;
		temp = z[j4 - 2];
		z[j4 - 2] = z[ipn4 - j4 - 2];
		z[ipn4 - j4 - 2] = temp;
		temp = z[j4 - 1];
		z[j4 - 1] = z[ipn4 - j4 - 5];
		z[ipn4 - j4 - 5] = temp;
		temp = z[j4];
		z[j4] = z[ipn4 - j4 - 4];
		z[ipn4 - j4 - 4] = temp;

	    }
	    if (n0 - i0 <= 4) {
		z[(n0 * 4) + pp - 1] = z[(i0 * 4) + pp - 1];
		z[(n0 * 4) - pp] = z[(i0 * 4) - pp];
	    }
	    mtemp1 = dmin2, mtemp2 = z[(n0 * 4) + pp - 1];
	    dmin2 = min(mtemp1, mtemp2);
	    mtemp1 = z[(n0 * 4) + pp - 1], mtemp2 = z[(i0 * 4) + pp - 1];
	    mtemp3 = min(mtemp1, mtemp2), mtemp2 = z[(i0 * 4) + pp + 3];
	    z[(n0 * 4) + pp - 1] = min(mtemp3, mtemp2);

	    mtemp1 = z[(n0 * 4) - pp], mtemp2 = z[(i0 * 4) - pp];
	    mtemp3 = min(mtemp1, mtemp2), mtemp2 = z[(i0 * 4) - pp + 4];
	    z[(n0 * 4) - pp] = min(mtemp3, mtemp2);
	    mtemp1 = qmax, mtemp2 = z[(i0 * 4) + pp - 3];
	    mtemp3 = max(mtemp1, mtemp2), mtemp2 = z[(i0 * 4) + pp + 1];
	    qmax = max(mtemp3, mtemp2);
	    *dmin = -Zero;
	}
    }
    mtemp1 = z[(n0 * 4) + pp - 1], mtemp2 = z[(n0 * 4) + pp - 9];
    mtemp3 = min(mtemp1, mtemp2), mtemp2 = dmin2 + z[(n0 * 4) - pp];
    if (*dmin < Zero || safmin * qmax < min(mtemp3, mtemp2)) {
//Choose a shift.
	Rlasq4(i0, n0, &z[0], pp, n0in, *dmin, dmin1, dmin2, dn, dn1, dn2, &tau, &ttype);
//Call dqds until DMIN > Zero
      L80:
	Rlasq5(i0, n0, &z[0], pp, tau, dmin, &dmin1, &dmin2, &dn, &dn1, &dn2, ieee);
	ndiv += n0 - i0 + 2;
	++(*iter);
//Check status.
	if (*dmin >= Zero && dmin1 > Zero) {
//Success.
	    goto L100;
	} else if (*dmin < Zero && dmin1 > Zero && z[(n0 - 1 * 4) - pp] < tol * (*sigma + dn1)
		   && abs(dn) < tol * *sigma) {
//Convergence hidden by negative DN.
	    z[(n0 - 1 * 4) - pp + 2] = Zero;
	    *dmin = Zero;
	    goto L100;
	} else if (*dmin < Zero) {
//TAU too big. Select new TAU and try again.
	    ++(nfail);
	    if (ttype < -22) {
//Failed twice. Play it safe.
		tau = Zero;
	    } else if (dmin1 > Zero) {
//Late failure. Gives excellent shift.
		tau = (tau + *dmin) * (One - eps * Two);
		ttype += -11;
	    } else {
//Early failure. Divide by Four
		tau *= Qurtr;
		ttype += -12;
	    }
	    goto L80;
	} else if (*dmin != *dmin) {
//NaN.
	    tau = Zero;
	    goto L80;
	} else {
//Possible underflow. Play it safe.
	    goto L90;
	}
    }
//Risk of underflow.
  L90:
    Rlasq6(i0, n0, &z[1], pp, dmin, &dmin1, &dmin2, &dn, &dn1, &dn2);
    ndiv += n0 - i0 + 2;
    ++(*iter);
    tau = Zero;

  L100:
    if (tau < *sigma) {
	*desig += tau;
	t = *sigma + *desig;
	*desig -= t - *sigma;
    } else {
	t = *sigma + tau;
	*desig = *sigma - (t - tau) + *desig;
    }
    *sigma = t;
    return;
}
