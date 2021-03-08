/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed6.cpp,v 1.12 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlaed6(INTEGER kniter, INTEGER orgati, REAL rho, REAL * d, REAL * z, REAL * finit, REAL * tau, INTEGER * info)
{
    REAL a, b, c, f;
    INTEGER i;
    REAL fc, df, ddf, lbd, eta, ubd, eps, base;
    INTEGER iter;
    REAL temp, temp1, temp2, temp3, temp4;
    INTEGER scale;
    INTEGER niter;
    REAL small1, small2, sminv1, sminv2;
    REAL dscale[3], sclfac, zscale[3], erretm, sclinv = 0;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Three = 3.0, Four = 4.0, Third;
    REAL mtemp1, mtemp2;

    *info = 0;
    if (orgati) {
	lbd = d[2];
	ubd = d[3];
    } else {
	lbd = d[1];
	ubd = d[2];
    }
    if (*finit < Zero) {
	lbd = Zero;
    } else {
	ubd = Zero;
    }

    niter = 1;
    *tau = Zero;
    if (kniter == 2) {
	if (orgati) {
	    temp = (d[3] - d[2]) / Two;
	    c = rho + z[1] / (d[1] - d[2] - temp);
	    a = c * (d[2] + d[3]) + z[2] + z[3];
	    b = c * d[2] * d[3] + z[2] * d[3] + z[3] * d[2];
	} else {
	    temp = (d[1] - d[2]) / Two;
	    c = rho + z[3] / (d[3] - d[2] - temp);
	    a = c * (d[1] + d[2]) + z[1] + z[2];
	    b = c * d[1] * d[2] + z[1] * d[2] + z[2] * d[1];
	}
	mtemp1 = abs(a), mtemp2 = abs(b);
	mtemp1 = max(mtemp1, mtemp2), mtemp2 = abs(c);
	temp = max(mtemp1, mtemp2);
	a /= temp;
	b /= temp;
	c /= temp;
	if (c == Zero) {
	    *tau = b / a;
	} else if (a <= Zero) {
	    *tau = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
	} else {
	    *tau = b * Two / (a + sqrt(abs(a * a - b * Four * c)));
	}
	if (*tau < lbd || *tau > ubd) {
	    *tau = (lbd + ubd) / Two;
	}
	if (d[1] == *tau || d[2] == *tau || d[3] == *tau) {
	    *tau = Zero;
	} else {
	    temp = *finit + *tau * z[1] / (d[1] * (d[1] - *tau)) + *tau * z[2] / (d[2] * (d[2] - *tau)) + *tau * z[3] / (d[3] * (d[3] - *tau));
	    if (temp <= Zero) {
		lbd = *tau;
	    } else {
		ubd = *tau;
	    }
	    if (abs(*finit) <= abs(temp)) {
		*tau = Zero;
	    }
	}
    }
//get machine parameters for possible scaling to avoid overflow
//modified by Sven: parameters SMALL1, SMINV1, SMALL2,
//SMINV2, EPS are not SAVEd anymore between one call to the
//others but recomputed at each call
    eps = Rlamch("Epsilon");
    base = Rlamch("Base");
    Third = One / Three;
    small1 = pow(Rlamch("SafMin"), Third);
    sminv1 = One / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;

//Determine if scaling of inputs necessary to avoid overflow
//when computing 1/TEMP**3

    if (orgati) {
	mtemp1 = abs(d[2] - *tau), mtemp2 = abs(d[3] - *tau);
	temp = min(mtemp1, mtemp2);
    } else {
	mtemp1 = abs(d[1] - *tau), mtemp2 = abs(d[2] - *tau);
	temp = min(mtemp1, mtemp2);
    }
    scale = MFALSE;
    if (temp <= small1) {
	scale = MTRUE;
	if (temp <= small2) {
//Scale up by power of radix nearest 1/SAFMIN**(2/3)
	    sclfac = sminv2;
	    sclinv = small2;
	} else {
//Scale up by power of radix nearest 1/SAFMIN**(1/3)
	    sclfac = sminv1;
	    sclinv = small1;
	}
//Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
	for (i = 0; i < 3; i++) {
	    dscale[i - 1] = d[i] * sclfac;
	    zscale[i - 1] = z[i] * sclfac;
	}
	*tau *= sclfac;
	lbd *= sclfac;
	ubd *= sclfac;
    } else {
//Copy D and Z to DSCALE and ZSCALE
	for (i = 0; i < 3; i++) {
	    dscale[i - 1] = d[i];
	    zscale[i - 1] = z[i];
	}
    }

    fc = Zero;
    df = Zero;
    ddf = Zero;
    for (i = 0; i < 3; i++) {
	temp = One / (dscale[i - 1] - *tau);
	temp1 = zscale[i - 1] * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / dscale[i - 1];
	df += temp2;
	ddf += temp3;
    }
    f = *finit + *tau * fc;

    if (abs(f) <= Zero) {
	goto L60;
    }
    if (f <= Zero) {
	lbd = *tau;
    } else {
	ubd = *tau;
    }
//Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
//                    scheme 
//It is not hard to see that
//       1) Iterations will go up monotonically
//              if FINIT < 0;
//       2) Iterations will go down monotonically
//              if FINIT > Zero
    iter = niter + 1;
    for (niter = iter; niter <= 40; ++niter) {
	if (orgati) {
	    temp1 = dscale[1] - *tau;
	    temp2 = dscale[2] - *tau;
	} else {
	    temp1 = dscale[0] - *tau;
	    temp2 = dscale[1] - *tau;
	}
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
	b = temp1 * temp2 * f;
	c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
	mtemp1 = abs(a), mtemp2 = abs(b);
	mtemp1 = max(mtemp1, mtemp2), mtemp2 = abs(c);
	temp = max(mtemp1, mtemp2);
	a /= temp;
	b /= temp;
	c /= temp;
	if (c == Zero) {
	    eta = b / a;
	} else if (a <= Zero) {
	    eta = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
	} else {
	    eta = b * Two / (a + sqrt(abs(a * a - b * Four * c)));
	}
	if (f * eta >= Zero) {
	    eta = -f / df;
	}
	*tau += eta;
	if (*tau < lbd || *tau > ubd) {
	    *tau = (lbd + ubd) / Two;
	}
	fc = Zero;
	erretm = Zero;
	df = Zero;
	ddf = Zero;
	for (i = 0; i < 3; i++) {
	    temp = One / (dscale[i - 1] - *tau);
	    temp1 = zscale[i - 1] * temp;
	    temp2 = temp1 * temp;
	    temp3 = temp2 * temp;
	    temp4 = temp1 / dscale[i - 1];
	    fc += temp4;
	    erretm += abs(temp4);
	    df += temp2;
	    ddf += temp3;
	}
	f = *finit + *tau * fc;
	erretm = (abs(*finit) + abs(*tau) * erretm) * 8. + abs(*tau) * df;
	if (abs(f) <= eps * erretm) {
	    goto L60;
	}
	if (f <= Zero) {
	    lbd = *tau;
	} else {
	    ubd = *tau;
	}

    }
    *info = 1;
  L60:

//Undo scaling
    if (scale) {
	*tau *= sclinv;
    }
    return;
}
