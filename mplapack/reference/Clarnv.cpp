/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarnv.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clarnv(INTEGER idist, INTEGER * iseed, INTEGER n, COMPLEX * x)
{
    INTEGER i;
    REAL u[128];
    INTEGER il, iv;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL TWOPI = 6.28;		//fake
    REAL tmp;
    COMPLEX ztmp;

    for (iv = 1; iv <= n; iv += 64) {
	il = min((INTEGER) 64, n - iv + 1);
//Call DLARUV to generate 2*IL real numbers from a uniform (0,1)
//distribution (2*IL <= LV)
	Rlaruv(&iseed[1], il * 2, u);
	if (idist == 1) {
//Copy generated numbers
	    for (i = 0; i < il; i++) {
		x[iv + i - 1] = u[(i * 2) - 2];
	    }
	} else if (idist == 2) {
//Convert generated numbers to uniform (-1,1) distribution
	    for (i = 0; i < il; i++) {
		x[iv + i - 1] = (u[(i * 2) - 2] * Two - One);
	    }
	} else if (idist == 3) {
//Convert generated numbers to normal (0,1) distribution
	    for (i = 0; i < il; i++) {
            tmp = sqrt(log(u[(i * 2) - 2]) * -Two);
		    x[iv + i - 1] = tmp * exp(ztmp);
	    }
	} else if (idist == 4) {
//Convert generated numbers to complex numbers uniformly
//distributed on the unit disk
	    for (i = 0; i < il; i++) {
		ztmp.real(Zero);
		ztmp.imag(TWOPI * u[(i * 2)]);
        tmp = sqrt(u[(i * 2) - 2]);
		x[iv + i - 1] = tmp * exp(ztmp);
	    }
	} else if (idist == 5) {
//Convert generated numbers to complex numbers uniformly
//distributed on the unit circle
	    for (i = 0; i < il; i++) {
		ztmp.real(Zero);
		ztmp.imag(TWOPI * u[(i * 2)]);
		x[iv + i - 1] = exp(ztmp);
	    }
	}
    }
    return;
}
