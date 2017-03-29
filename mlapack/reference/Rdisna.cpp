/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rdisna.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE  1
#define MFALSE 0

void Rdisna(const char *job, INTEGER m, INTEGER n, REAL * d, REAL * sep, INTEGER * info)
{

    INTEGER i, k;
    REAL eps;
    INTEGER decr, left, incr, sing, eigen;
    REAL anorm;
    INTEGER right;
    REAL oldgap, safmin;
    REAL newgap, thresh;
    REAL Zero = 0.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    eigen = Mlsame(job, "E");
    left = Mlsame(job, "L");
    right = Mlsame(job, "R");
    sing = left || right;
    if (eigen) {
	k = m;
    } else if (sing) {
	k = min(m, n);
    }
    if (!eigen && !sing) {
	*info = -1;
    } else if (m < 0) {
	*info = -2;
    } else if (k < 0) {
	*info = -3;
    } else {
	incr = MTRUE;
	decr = MTRUE;
	for (i = 0; i < k - 1; i++) {
	    if (incr) {
		incr = incr && d[i] <= d[i + 1];
	    }
	    if (decr) {
		decr = decr && d[i] >= d[i + 1];
	    }

	}
	if (sing && k > 0) {
	    if (incr) {
		incr = incr && Zero <= d[1];
	    }
	    if (decr) {
		decr = decr && d[k] >= Zero;
	    }
	}
	if (!(incr || decr)) {
	    *info = -4;
	}
    }
    if (*info != 0) {
	Mxerbla("Rdisna", -(*info));
	return;
    }
//Quick return if possible
    if (k == 0)
	return;

//Compute reciprocal condition numbers
    if (k == 1) {
	sep[1] = Rlamch("O");
    } else {
	oldgap = (abs(d[2] - d[1]));
	sep[1] = oldgap;
	for (i = 1; i < k - 1; i++) {
	    newgap = abs(d[i + 1] - d[i]);
	    sep[i] = min(oldgap, newgap);
	    oldgap = newgap;

	}
	sep[k] = oldgap;
    }
    if (sing) {
	if ((left && m > n) || (right && m < n)) {
	    if (incr) {
		sep[1] = min(sep[1], d[1]);
	    }
	    if (decr) {
		mtemp1 = sep[k], mtemp2 = d[k];
		sep[k] = min(mtemp1, mtemp2);
	    }
	}
    }
//Ensure that reciprocal condition numbers are not less than
//threshold, in order to limit the size of the error bound
    eps = Rlamch("E");
    safmin = Rlamch("S");
    mtemp1 = abs(d[1]);
    mtemp2 = abs(d[k]);
    anorm = max(mtemp1, mtemp2);
    if (anorm == Zero) {
	thresh = eps;
    } else {
	mtemp1 = eps * anorm;
	thresh = max(mtemp1, safmin);
    }
    for (i = 0; i < k; i++) {
	mtemp1 = sep[i];
	sep[i] = max(mtemp1, thresh);
    }
    return;
}
