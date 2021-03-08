/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rpttrf.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rpttrf(INTEGER n, REAL * d, REAL * e, INTEGER * info)
{
    REAL Zero = 0.0;
    INTEGER i, i4 = 0;
    REAL ei;

    *info = 0;
    if (n < 0) {
	*info = -1;
	Mxerbla("Rpttrf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

//Compute the L*D*L' (or U'*D*U) factorization of A.
    for (i = 0; i < (n - 1) % 4; i++) {
	if (d[i] <= Zero) {
	    *info = i;
	    return;
	}
	ei = e[i];
	e[i] = ei / d[i];
	d[i + 1] -= e[i] * ei;
    }

    for (i = i4 + 1; i < n - 4; i += 4) {
//Drop out of the loop if d(i) <= 0: the matrix is not positive
//definite.

	if (d[i] <= Zero) {
	    *info = i;
	    return;
	}
//Solve for e(i) and d(i+1).
	ei = e[i];
	e[i] = ei / d[i];
	d[i + 1] -= e[i] * ei;

	if (d[i + 1] <= Zero) {
	    *info = i + 1;
	    return;
	}
//Solve for e(i+1) and d(i+2).
	ei = e[i + 1];
	e[i + 1] = ei / d[i + 1];
	d[i + 2] -= e[i + 1] * ei;

	if (d[i + 2] <= Zero) {
	    *info = i + 2;
	    return;
	}
//Solve for e(i+2) and d(i+3).

	ei = e[i + 2];
	e[i + 2] = ei / d[i + 2];
	d[i + 3] -= e[i + 2] * ei;

	if (d[i + 3] <= Zero) {
	    *info = i + 3;
	    return;
	}
//Solve for e(i+3) and d(i+4).

	ei = e[i + 3];
	e[i + 3] = ei / d[i + 3];
	d[i + 4] -= e[i + 3] * ei;
    }

//Check d(n) for positive definiteness.

    if (d[n] <= Zero) {
	*info = n;
    }
    return;
}
