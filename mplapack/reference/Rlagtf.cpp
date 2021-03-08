/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlagtf.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlagtf(INTEGER n, REAL * a, REAL lambda, REAL * b, REAL * c, REAL * tol, REAL * d, INTEGER * in, INTEGER * info)
{
    INTEGER k;
    REAL tl, eps, piv1, piv2, temp, mult, scale1, scale2;
    REAL Zero = 0.0;

    *info = 0;
    if (n < 0) {
	*info = -1;
	Mxerbla("Rlagtf", -(*info));
	return;
    }

    if (n == 0)
	return;

    a[0] -= lambda;
    in[n] = 0;
    if (n == 1) {
	if (a[0] == Zero) {
	    in[1] = 1;
	}
	return;
    }

    eps = Rlamch("Epsilon");

    tl = max(*tol, eps);
    scale1 = abs(a[1]) + abs(b[1]);
    for (k = 0; k < n - 1; k++) {
	a[k + 1] -= lambda;
	scale2 = abs(c[k]) + abs(a[k + 1]);
	if (k < n - 1) {
	    scale2 += abs(b[k + 1]);
	}
	if (a[k] == Zero) {
	    piv1 = Zero;
	} else {
	    piv1 = abs(a[k]) / scale1;
	}
	if (c[k] == Zero) {
	    in[k] = 0;
	    piv2 = Zero;
	    scale1 = scale2;
	    if (k < n - 1) {
		d[k] = Zero;
	    }
	} else {
	    piv2 = abs(c[k]) / scale2;
	    if (piv2 <= piv1) {
		in[k] = 0;
		scale1 = scale2;
		c[k] /= a[k];
		a[k + 1] -= c[k] * b[k];
		if (k < n - 1) {
		    d[k] = Zero;
		}
	    } else {
		in[k] = 1;
		mult = a[k] / c[k];
		a[k] = c[k];
		temp = a[k + 1];
		a[k + 1] = b[k] - mult * temp;
		if (k < n - 1) {
		    d[k] = b[k + 1];
		    b[k + 1] = -mult * d[k];
		}
		b[k] = temp;
		c[k] = mult;
	    }
	}
	if (max(piv1, piv2) <= tl && in[n] == 0) {
	    in[n] = k;
	}

    }
    if (abs(a[n]) <= scale1 * tl && in[n] == 0) {
	in[n] = n;
    }
    return;
}
