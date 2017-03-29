/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claic1.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Claic1(INTEGER job, INTEGER j, COMPLEX * x, REAL sest, COMPLEX * w, COMPLEX gamma, REAL * sestpr, COMPLEX * s, COMPLEX * c)
{
    REAL b, t, s1, s2, scl, eps, tmp;
    COMPLEX sine;
    REAL test, zeta1, zeta2;
    COMPLEX alpha;
    REAL norma;
    REAL absgam, absalp;
    COMPLEX cosine;
    REAL absest;
    REAL Zero = 0.0, One = 1.0, Half = 0.5, Two = 2.0, Four = 4.0;
    REAL mtemp1, mtemp2;
    COMPLEX ctemp1;

    eps = Rlamch("Epsilon");
    alpha = Cdotc(j, &x[0], 1, &w[1], 1);

    absalp = abs(alpha);
    absgam = abs(gamma);
    absest = abs(sest);
    if (job == 1) {
//Estimating largest singular value
//special cases

	if (sest == Zero) {
	    s1 = max(absgam, absalp);
	    if (s1 == Zero) {
		*s = Zero;
		*c = One;
		*sestpr = Zero;
	    } else {
		*s = alpha / s1;
		*c = gamma / s1;
		mtemp1 = (*s * (conj(*s)) + *c * (conj(*c))).real();	//always real
		tmp = sqrt(mtemp1);
		*s = *s / tmp;
		*c = *c / tmp;
		*sestpr = s1 * tmp;
	    }
	    return;
	} else if (absgam <= eps * absest) {
	    *s = One;
	    *c = Zero;
	    tmp = max(absest, absalp);
	    s1 = absest / tmp;
	    s2 = absalp / tmp;
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
	    return;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		*s = One;
		*c = Zero;
		*sestpr = s2;
	    } else {
		*s = Zero;
		*c = One;
		*sestpr = s1;
	    }
	    return;
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		scl = sqrt(tmp * tmp + One);
		*sestpr = s2 * scl;
		*s = alpha / s2 / scl;
		*c = gamma / s2 / scl;
	    } else {
		tmp = s2 / s1;
		scl = sqrt(tmp * tmp + One);
		*sestpr = s1 * scl;
		*s = alpha / s1 / scl;
		*c = gamma / s1 / scl;
	    }
	    return;
	} else {
//normal case
	    zeta1 = absalp / absest;
	    zeta2 = absgam / absest;
	    b = (One - zeta1 * zeta1 - zeta2 * zeta2) * Half;
	    *c = zeta1 * zeta1;
	    if (b > Zero) {
		t = (*c / (b + sqrt(b * b + *c))).real();
	    } else {
		t = (sqrt(b * b + *c) - b).real();
	    }

	    sine = -(alpha / absest) / t;
	    cosine = -(gamma / absest) / (t + One);
	    tmp = sqrt(sine * conj(sine) + cosine * conj(cosine)).real();
	    *s = sine / tmp;
	    *c = cosine / tmp;
	    *sestpr = sqrt(t + One) * absest;
	    return;
	}

    } else if (job == 2) {
//Estimating smallest singular value
//special cases
	if (sest == Zero) {
	    *sestpr = Zero;
	    if (max(absgam, absalp) == Zero) {
		sine = One;
		cosine = Zero;
	    } else {
		sine = -conj(gamma);
		cosine = conj(alpha);
	    }
	    mtemp1 = abs(sine), mtemp2 = abs(cosine);
	    s1 = max(mtemp1, mtemp2);
	    *s = sine / s1;
	    *c = cosine / s1;
	    tmp = sqrt(*s * conj(*s) + *c * conj(*c)).real();
	    *s = *s / tmp;
	    *c = *c / tmp;
	    return;
	} else if (absgam <= eps * absest) {
	    *s = Zero;
	    *c = One;
	    *sestpr = absgam;
	    return;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		*s = Zero;
		*c = One;
		*sestpr = s1;
	    } else {
		*s = One;
		*c = Zero;
		*sestpr = s2;
	    }
	    return;
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		scl = sqrt(tmp * tmp + One);
		*sestpr = absest * (tmp / scl);
		*s = -(conj(gamma) / s2) / scl;
		*c = conj(alpha) / s2 / scl;
	    } else {
		tmp = s2 / s1;
		scl = sqrt(tmp * tmp + One);
		*sestpr = absest / scl;
		*s = -(conj(gamma) / s1) / scl;
		*c = conj(alpha) / s1 / scl;
	    }
	    return;
	} else {
//normal case
	    zeta1 = absalp / absest;
	    zeta2 = absgam / absest;
	    mtemp1 = zeta1 * zeta1 + One + zeta1 * zeta2;
	    mtemp2 = zeta1 * zeta2 + zeta2 * zeta2;
	    norma = max(mtemp1, mtemp2);
//See if root is closer to zero or to ONE
	    test = (zeta1 - zeta2) * Two * (zeta1 + zeta2) + One;
	    if (test >= Zero) {
//root is close to zero, compute directly
		b = (zeta1 * zeta1 + zeta2 * zeta2 + One) * Half;
		*c = zeta2 * zeta2;
		t = (*c / (b + sqrt(abs(b * b - *c)))).real();
		sine = alpha / absest / (One - t);
		cosine = -(gamma / absest) / t;
		*sestpr = sqrt(t + eps * Four * eps * norma) * absest;
	    } else {
//root is closer to ONE, shift by that amount
		b = (zeta2 * zeta2 + zeta1 * zeta1 - One) * Half;
		*c = zeta1 * zeta1;
		if (b >= Zero) {
		    t = (-(*c) / (b + sqrt(b * b + *c))).real();
		} else {
		    t = (b - sqrt(b * b + *c)).real();
		}
		sine = -(alpha / absest) / t;
		cosine = -(gamma / absest) / (t + One);
		*sestpr = sqrt(t + One + eps * Four * eps * norma) * absest;
	    }
	    tmp = sqrt(sine * conj(sine) + cosine * conj(cosine)).real();
	    *s = sine / tmp;
	    *c = cosine / tmp;
	    return;
	}
    }
    return;
}
