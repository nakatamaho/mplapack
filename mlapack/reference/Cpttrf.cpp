/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cpttrf.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cpttrf(INTEGER n, REAL * d, COMPLEX * e, INTEGER * info)
{
    REAL f, g;
    INTEGER i, i4;
    REAL eii, eir;
    REAL Zero = 0.0;

//Test the input parameters.
    *info = 0;
    if (n < 0) {
	*info = -1;
	Mxerbla("Cpttrf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Compute the L*D*L' (or U'*D*U) factorization of A.
    i4 = (n - 1) % 4;
    for (i = 0; i < i4; i++) {
	if (d[i] <= Zero) {
	    *info = i;
	    goto L30;
	}
	eir = e[i].real();
	eii = e[i].imag();
	f = eir / d[i];
	g = eii / d[i];
	e[i] = f;
	d[i + 1] = d[i + 1] - f * eir - g * eii;
    }
    for (i = i4 + 1; i <= n - 4; i = i + 4) {
//Drop out of the loop if d(i) <= 0: the matrix is not positive
//definite.
	if (d[i] <= Zero) {
	    *info = i;
	    goto L30;
	}
//Solve for e(i) and d(i+1).
	eir = e[i].real();
	eii = e[i].imag();
	f = eir / d[i];
	g = eii / d[i];
	e[i] = f;
	d[i + 1] = d[i + 1] - f * eir - g * eii;
	if (d[i + 1] <= Zero) {
	    *info = i + 1;
	    goto L30;
	}
//Solve for e(i+1) and d(i+2).
	eir = e[i + 1].real();
	eii = e[i + 1].imag();
	f = eir / d[i + 1];
	g = eii / d[i + 1];
	e[i + 1] = f;
	d[i + 2] = d[i + 2] - f * eir - g * eii;
	if (d[i + 2] <= Zero) {
	    *info = i + 2;
	    goto L30;
	}
//Solve for e(i+2) and d(i+3).
	eir = e[i + 2].real();
	eii = e[i + 2].imag();
	f = eir / d[i + 2];
	g = eii / d[i + 2];
	e[i + 2] = f;
	d[i + 3] = d[i + 3] - f * eir - g * eii;
	if (d[i + 3] <= Zero) {
	    *info = i + 3;
	    goto L30;
	}
//Solve for e(i+3) and d(i+4).
	eir = e[i + 3].real();
	eii = e[i + 3].imag();
	f = eir / d[i + 3];
	g = eii / d[i + 3];
	e[i + 3] = f;
	d[i + 4] = d[i + 4] - f * eir - g * eii;
    }
//Check d(n) for positive definiteness.
    if (d[n] <= Zero) {
	*info = n;
    }
  L30:
    return;
}
