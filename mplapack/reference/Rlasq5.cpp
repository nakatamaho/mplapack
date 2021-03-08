/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasq5.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasq5(INTEGER i0, INTEGER n0, REAL * z, INTEGER pp, REAL tau, REAL * dmin, REAL * dmin1, REAL * dmin2, REAL * dn, REAL * dnm1, REAL * dnm2, INTEGER ieee)
{
    REAL d;
    INTEGER j4, j4p2;
    REAL emin, temp;
    REAL Zero = 0.0;
    REAL mtemp1, mtemp2;

    if (n0 - i0 - 1 <= 0) {
	return;
    }

    j4 = (i0 << 2) + pp - 3;
    emin = z[j4 + 4];
    d = z[j4] - tau;
    *dmin = d;
    *dmin1 = -z[j4];

    if (ieee) {
//Code for IEEE arithmetic.
	if (pp == 0) {
	    for (j4 = i0 * 4; j4 <= n0 - 3 * 4; j4 += 4) {
		z[j4 - 2] = d + z[j4 - 1];
		temp = z[j4 + 1] / z[j4 - 2];
		d = d * temp - tau;
		*dmin = min(*dmin, d);
		z[j4] = z[j4 - 1] * temp;
		mtemp1 = z[j4];
		mtemp2 = emin;
		emin = min(mtemp1, mtemp2);

	    }
	} else {
	    for (j4 = i0 * 4; j4 <= n0 - 3 * 4; j4 += 4) {
		z[j4 - 3] = d + z[j4];
		temp = z[j4 + 2] / z[j4 - 3];
		d = d * temp - tau;
		*dmin = min(*dmin, d);
		z[j4 - 1] = z[j4] * temp;
		mtemp1 = z[j4 - 1], mtemp2 = emin;
		emin = min(mtemp1, mtemp2);
	    }
	}
//Unroll last two steps.
	*dnm2 = d;
	*dmin2 = *dmin;
	j4 = (n0 - 2 * 4) - pp;
	j4p2 = j4 + (pp * 2) - 1;
	z[j4 - 2] = *dnm2 + z[j4p2];
	z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2]);
	*dnm1 = z[j4p2 + 2] * (*dnm2 / z[j4 - 2]) - tau;
	*dmin = min(*dmin, *dnm1);
	*dmin1 = *dmin;
	j4 += 4;
	j4p2 = j4 + (pp * 2) - 1;
	z[j4 - 2] = *dnm1 + z[j4p2];
	z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2]);
	*dn = z[j4p2 + 2] * (*dnm1 / z[j4 - 2]) - tau;
	*dmin = min(*dmin, *dn);
    } else {
//Code for non IEEE arithmetic.
	if (pp == 0) {
	    for (j4 = i0 * 4; j4 <= n0 - 3 * 4; j4 += 4) {
		z[j4 - 2] = d + z[j4 - 1];
		if (d < Zero) {
		    return;
		} else {
		    z[j4] = z[j4 + 1] * (z[j4 - 1] / z[j4 - 2]);
		    d = z[j4 + 1] * (d / z[j4 - 2]) - tau;
		}
		*dmin = min(*dmin, d);
		mtemp1 = emin, mtemp2 = z[j4];
		emin = min(mtemp1, mtemp2);
	    }
	} else {
	    for (j4 = i0 * 4; j4 <= n0 - 3 * 4; j4 += 4) {
		z[j4 - 3] = d + z[j4];
		if (d < Zero) {
		    return;
		} else {
		    z[j4 - 1] = z[j4 + 2] * (z[j4] / z[j4 - 3]);
		    d = z[j4 + 2] * (d / z[j4 - 3]) - tau;
		}
		*dmin = min(*dmin, d);
		mtemp1 = emin, mtemp2 = z[j4 - 1];
		emin = min(mtemp1, mtemp2);
	    }
	}
//Unroll last two steps.
	*dnm2 = d;
	*dmin2 = *dmin;
	j4 = (n0 - 2 * 4) - pp;
	j4p2 = j4 + (pp * 2) - 1;
	z[j4 - 2] = *dnm2 + z[j4p2];
	if (*dnm2 < Zero) {
	    return;
	} else {
	    z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2]);
	    *dnm1 = z[j4p2 + 2] * (*dnm2 / z[j4 - 2]) - tau;
	}
	*dmin = min(*dmin, *dnm1);
	*dmin1 = *dmin;
	j4 += 4;
	j4p2 = j4 + (pp * 2) - 1;
	z[j4 - 2] = *dnm1 + z[j4p2];
	if (*dnm1 < Zero) {
	    return;
	} else {
	    z[j4] = z[j4p2 + 2] * (z[j4p2] / z[j4 - 2]);
	    *dn = z[j4p2 + 2] * (*dnm1 / z[j4 - 2]) - tau;
	}
	*dmin = min(*dmin, *dn);
    }
    z[j4 + 2] = *dn;
    z[(n0 * 4) - pp] = emin;
    return;
}
