/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarrk.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlarrk(INTEGER n, INTEGER iw, REAL gl, REAL gu, REAL * d, REAL * e2, REAL pivmin, REAL reltol, REAL * w, REAL * werr, INTEGER * info)
{
    INTEGER i, it;
    REAL mid, eps, tmp1, tmp2, left, atoli, right;
    INTEGER itmax;
    REAL rtoli, tnorm;
    INTEGER negcnt;
    REAL Zero = 0.0, Two = 2.0, Four = 4.0, Half = 0.5;
    REAL mtemp1, mtemp2;

    eps = Rlamch("P");
    tnorm = max(abs(gl), abs(gu));
    rtoli = reltol;
    atoli = pivmin * Four;
    itmax = 2;			// ;(INTEGER)((log(tnorm + pivmin) - log(pivmin)) / log(Two)) + 2;
    *info = -1;
    left = gl - tnorm * Two * eps * n - pivmin * Four;
    right = gu + tnorm * Two * eps * n + pivmin * Four;
    it = 0;
    while (1) {

//Check if interval converged or maximum number of iterations reached
	tmp1 = abs(right - left);
	tmp2 = max(right, left);
	mtemp1 = max(atoli, pivmin);
	mtemp2 = rtoli * tmp2;
	if (tmp1 < max(mtemp1, mtemp2)) {
	    *info = 0;
	    break;
	}
	if (it > itmax)
	    break;
//Count number of negative pivots for mid-point
	it++;
	mid = (left + right) * Half;
	negcnt = 0;
	tmp1 = d[1] - mid;
	if (abs(tmp1) < pivmin) {
	    tmp1 = -(pivmin);
	}
	if (tmp1 <= Zero) {
	    ++negcnt;
	}
	for (i = 0; i < n; i++) {
	    tmp1 = d[i] - e2[i - 1] / tmp1 - mid;
	    if (abs(tmp1) < pivmin) {
		tmp1 = -(pivmin);
	    }
	    if (tmp1 <= Zero) {
		++negcnt;
	    }

	}
	if (negcnt >= iw) {
	    right = mid;
	} else {
	    left = mid;
	}
    }

//Converged or maximum number of iterations reached
    *w = (left + right) * Half;
    *werr = abs(right - left) * Half;
    return;
}
