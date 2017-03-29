/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlanv2.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlanv2(REAL * a, REAL * b, REAL * c, REAL * d, REAL * rt1r, REAL * rt1i, REAL * rt2r, REAL * rt2i, REAL * cs, REAL * sn)
{

    REAL p, z, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, scale, bcmax, bcmis, sigma;

    REAL One = 1.0, Half = 0.5, Zero = 0.0, Four = 4.0;
    REAL mtemp1, mtemp2;

//Parameters of the rotation matrix.
    eps = Rlamch("P");
    if (*c == Zero) {
	*cs = One;
	*sn = Zero;
	goto L10;
    } else if (*b == Zero) {
//Swap rows and columns
	*cs = Zero;
	*sn = One;
	temp = *d;
	*d = *a;
	*a = temp;
	*b = -(*c);
	*c = Zero;
	goto L10;
    } else if (*a - *d == Zero && sign(One, *b) != sign(One, *c)) {
	*cs = One;
	*sn = Zero;
	goto L10;
    } else {

	temp = *a - *d;
	p = temp * Half;
	bcmax = max(*b, *c);
	bcmis = min(*b, *c) * sign(One, *b) * sign(One, *c);
	mtemp1 = abs(p);
	mtemp2 = bcmax;
	scale = max(mtemp1, mtemp2);
	z = p / scale * p + bcmax / scale * bcmis;

//If Z is of the order of the machine accuracy, postpone the
//decision on the nature of eigenvalues
	if (z >= eps * Four) {

//Real eigenvalues. Compute A and D.
	    mtemp1 = sqrt(scale) * sqrt(z);
	    z = p + sign(mtemp1, p);
	    *a = *d + z;
	    *d -= bcmax / z * bcmis;

//Compute B and the rotation matrix
	    tau = Rlapy2(*c, z);
	    *cs = z / tau;
	    *sn = *c / tau;
	    *b -= *c;
	    *c = Zero;
	} else {

//Complex eigenvalues, or real (almost) equal eigenvalues.
//Make diagonal elements equal.
	    sigma = *b + *c;
	    tau = Rlapy2(sigma, temp);
	    *cs = sqrt((abs(sigma) / tau + One) * Half);
	    *sn = -(p / (tau * *cs)) * sign(One, sigma);

//Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
//        [ CC  DD ]   [ C  D ] [ SN  CS ]

	    aa = *a * *cs + *b * *sn;
	    bb = -(*a) * *sn + *b * *cs;
	    cc = *c * *cs + *d * *sn;
	    dd = -(*c) * *sn + *d * *cs;

//Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
//        [ C  D ]   [-SN  CS ] [ CC  DD ]

	    *a = aa * *cs + cc * *sn;
	    *b = bb * *cs + dd * *sn;
	    *c = -aa * *sn + cc * *cs;
	    *d = -bb * *sn + dd * *cs;

	    temp = (*a + *d) * Half;
	    *a = temp;
	    *d = temp;

	    if (*c != Zero) {
		if (*b != Zero) {
		    if (sign(One, *b) == sign(One, *c)) {
//Real eigenvalues: reduce to upper triangular form
			sab = sqrt((abs(*b)));
			sac = sqrt((abs(*c)));
			mtemp1 = sab * sac;
			p = sign(mtemp1, *c);
			tau = One / sqrt(abs(*b + *c));
			*a = temp + p;
			*d = temp - p;
			*b -= *c;
			*c = Zero;
			cs1 = sab * tau;
			sn1 = sac * tau;
			temp = *cs * cs1 - *sn * sn1;
			*sn = *cs * sn1 + *sn * cs1;
			*cs = temp;
		    }
		} else {
		    *b = -(*c);
		    *c = Zero;
		    temp = *cs;
		    *cs = -(*sn);
		    *sn = temp;
		}
	    }
	}

    }

  L10:

//Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
    *rt1r = *a;
    *rt2r = *d;
    if (*c == Zero) {
	*rt1i = Zero;
	*rt2i = Zero;
    } else {
	*rt1i = sqrt((abs(*b))) * sqrt((abs(*c)));
	*rt2i = -(*rt1i);
    }
    return;
}
