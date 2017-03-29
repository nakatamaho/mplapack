/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claesy.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Claesy(COMPLEX a, COMPLEX b, COMPLEX c, COMPLEX * rt1, COMPLEX * rt2, COMPLEX * evscal, COMPLEX * cs1, COMPLEX * sn1)
{
    COMPLEX s, t;
    REAL z;
    COMPLEX tmp;
    REAL babs, tabs, evnorm;
    REAL Zero = 0.0, Half = 0.5, One = 1.0;
    REAL Thresh = 0.1;

//Special case:  The matrix is actually diagonal.
//To avoid divide by zero later, we treat this case separately.
    if (abs(b) == Zero) {
	*rt1 = a;
	*rt2 = c;
	if (abs(*rt1) < abs(*rt2)) {
	    tmp = *rt1;
	    *rt1 = *rt2;
	    *rt2 = tmp;
	    *cs1 = Zero;
	    *sn1 = One;
	} else {
	    *cs1 = One;
	    *sn1 = Zero;
	}
    } else {
//Compute the eigenvalues and eigenvectors.
//The characteristic equation is
//   lambda **2 - (A+C) lambda + (A*C - B*B)
//and we solve it using the quadratic formula.
	s = (a + c) * Half;
	t = (a - c) * Half;
//Take the square root carefully to avoid over/under flow.
	babs = abs(b);
	tabs = abs(t);
	z = max(babs, tabs);
	if (z > Zero) {
	    t = z * sqrt((t / z) * (t / z) + (b / z) * (b / z));
	}
//Compute the two eigenvalues.  RT1 and RT2 are exchanged
//if necessary so that RT1 will have the greater magnitude.
	*rt1 = s + t;
	*rt2 = s - t;
	if (abs(*rt1) < abs(*rt2)) {
	    tmp = *rt1;
	    *rt1 = *rt2;
	    *rt2 = tmp;
	}
//Choose CS1 = 1 and SN1 to satisfy the first equation, then
//scale the components of this eigenvector so that the matrix
//of eigenvectors X satisfies  X * X' = I .  (No scaling is
//done if the norm of the eigenvalue matrix is less than THRESH.)
	*sn1 = (*rt1 - a) / b;
	tabs = abs(*sn1);
	if (tabs > One) {
	    t = tabs * sqrt((One / tabs) * (One / tabs) + (*sn1 / tabs) * (*sn1 / tabs));
	} else {
	    t = sqrt(One + *sn1 * *sn1);
	}
	evnorm = abs(t);
	if (evnorm >= Thresh) {
	    *evscal = One / t;
	    *cs1 = *evscal;
	    *sn1 = *sn1 * *evscal;
	} else {
	    *evscal = Zero;
	}
    }
    return;
}
