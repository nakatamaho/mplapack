/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaev2.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

// http://www.netlib.org/lapack/double/dlaev2.f

#include <mpblas.h>
#include <mplapack.h>
#include <stdio.h>		//for printf. shall be removed

void Rlaev2(REAL a, REAL b, REAL c, REAL * rt1, REAL * rt2, REAL * cs1, REAL * sn1)
{
    REAL ab, acmn, acmx, acs, adf;
    REAL cs, ct, df, rt, sm, tb, tn;
    REAL zero, one, two, half;
    INTEGER sgn1, sgn2;

    zero = 0.0;
    one = 1.0;
    two = 2.0;
    half = 0.5;

    sm = a + c;
    df = a - c;
    adf = abs(df);
    tb = b + b;
    ab = abs(tb);

    if (abs(a) > abs(c)) {
	acmx = a;
	acmn = c;
    } else {
	acmx = c;
	acmn = a;
    }
    if (adf > ab) {
	rt = adf * sqrt(one + (ab / adf) * (ab / adf));
    } else if (adf < ab) {
	rt = ab * sqrt(one + (adf / ab) * (adf / ab));
    } else {
//Includes case AB=ADF=0
	rt = ab * sqrt(two);
    }
    if (sm < zero) {
	*rt1 = half * (sm - rt);
	sgn1 = -1;
//Order of execution important.
//To get fully accurate smaller eigenvalue,
//next line needs to be executed in higher precision.
	*rt2 = (acmx / (*rt1)) * acmn - (b / (*rt1)) * b;
    } else if (sm > zero) {
	*rt1 = half * (sm + rt);
	sgn1 = 1;
//Order of execution important.
//To get fully accurate smaller eigenvalue,
//next line needs to be executed in higher precision.
	*rt2 = (acmx / (*rt1)) * acmn - (b / (*rt1)) * b;
    } else {
//Includes case RT1 = RT2 = 0
	*rt1 = half * rt;
	*rt2 = -1.0 * half * rt;
	sgn1 = 1;
    }
//Compute the eigenvector
    if (df >= zero) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = one / sqrt(one + ct * ct);
	*cs1 = ct * (*sn1);
    } else {
	if (ab == zero) {
	    *cs1 = one;
	    *sn1 = zero;
	} else {
	    tn = -cs / tb;
	    *cs1 = one / sqrt(one + tn * tn);
	    *sn1 = tn * (*cs1);
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return;
}
